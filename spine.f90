module spine
    implicit none
    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer, parameter :: eout = 30, cvvout = 31
    integer :: Natom, Nbath, naccumulated, window_width, ncvvout
    real*8 :: rcmin, tstart, tstop, dtout, dt, bl, bl2,rho, kT
    real*8, dimension(:), allocatable :: y, Qbath, xi, vxi
    real*8, dimension(:,:), allocatable :: r0,  r, p, v, rshift, r0equil,  Cvvcur, Cvvold
    real*8, dimension(:,:,:), allocatable :: Qnk, Meff, invMeff, v0tau, v0s, p0
    real*8 :: lastepot
    character(256) :: stem

    interface
        subroutine nose_hoover_chain(p, Ekin, kT, xi, vxi, Q, dt, ne)
            real*8, intent(inout) :: p(:,:), Ekin, xi(:), vxi(:)
            real*8, intent(in) :: kT, Q(:), dt
            integer, intent(in) :: ne
        end subroutine nose_hoover_chain
    end interface
contains

function matsqrt(M)
    implicit none
    real*8, intent(in) :: M(:,:)
    real*8 :: matsqrt(size(M,1),size(M,2))
    real*8 :: U(size(M,1),size(M,1)), U1(size(M,1),size(M,1)), W(size(M,1)), FV1(size(M,1)), FV2(size(M,1))
    integer :: i, j, ierr, N

    N = size(M,1)
    call RS(N,N,M,W,1,U,FV1,FV2,ierr)
    U1 = transpose(U)
    do j=1,N
        U1(:,j) = sqrt(W)*U1(:,j)
    end do
    matsqrt = matmul(U, U1)
end function

subroutine initial_momenta(kT, Meff, p)
    use utils
    implicit none
    real*8, intent(in) :: kT, Meff(:,:,:)
    real*8, intent(out) :: p(:, :)
    real*8 :: U(3,3), W(3), FV1(3), FV2(3), p3(3)
    integer :: i, j, ierr
    
    do i=1,size(Meff, 3)
        call RS(3,3,Meff(:,:,i),W,1,U,FV1,FV2,ierr)
        do j=1,3
            p3(j) = gaussran(sqrt(W(j)*kT), 0.0d0)
        end do
        p(:,i) = matmul(U, p3)
    end do
end subroutine


subroutine verletstep(r, p, f, epot, dt)
    use vgw
    implicit none
    real*8, intent(in) :: dt
    real*8, intent(inout) ::r(:,:), p(:,:), f(:,:)
    real*8, intent(out) :: epot
    real*8 :: Ueff, v3(3), p3(3), Ekin,fcm(3)
    integer :: i, Nq, Nqp

    p = p + 0.5*dt*f
    do i=1,Natom
        r(:,i) = r(:,i) + dt*matmul(invMeff(:,:,i), p(:,i))
    end do
    call vgw1(r, Ueff, 1.0/kT, 0.0d0, y, Meff, invMeff)
    f = reshape( 2.0*kT*y(2+18*Natom:1+21*Natom), (/ 3, Natom /))
    fcm = sum(f,2)/Natom
    do i=1,Natom
        f(:,i) = f(:,i) - fcm
    end do
    p = p + 0.5*dt*f

    epot = -2.0*kT*y(1)
end subroutine


subroutine total_ekin(ekin)
    implicit none
    real*8, intent(out) :: ekin
    real*8 :: v3(3), p3(3)
    integer :: j

    ekin = 0
    do j=1,Natom
        !p3 = qp(1+3*(Natom+j-1):3*(Natom+j))
        v3 = matmul(invMeff(:,:,j), p(:,j))
        ekin = ekin + dot_product(p(:,j), v3)
    end do
    ekin = 0.5*ekin
end subroutine

subroutine velocity_autocorrelation()
    use vgw
    implicit none
    integer :: i,j
    real * 8 :: vs(3,Natom)

    do j=1,Natom
        v(:,j) = matmul(invMeff(:,:,j), p(:,j))
        vs(:,j) = matmul(matsqrt(invMeff(:,:,j)), p(:,j))/sqrt(mass)
    end do

    naccumulated = naccumulated + 1
    call unpack_Qnk(y, Qnk)
    do j=1,Natom
        v0tau(:,j,naccumulated) = matmul(Qnk(:,:,j), v(:,j))
    end do
    v0s(:,:,naccumulated) = vs
    p0(:,:,naccumulated) = p/mass

    do i=1,naccumulated
        Cvvcur(1,naccumulated-i+1) = Cvvcur(1,naccumulated-i+1) + sum(v0tau(:,:,i)*v)
        Cvvcur(2,naccumulated-i+1) = Cvvcur(2,naccumulated-i+1) + sum(v0s(:,:,i)*vs)
        Cvvcur(3,naccumulated-i+1) = Cvvcur(3,naccumulated-i+1) + sum(p0(:,:,i)*v)
    end do

    j = window_width
    do i=naccumulated+1,window_width
        Cvvold(1,j) =  Cvvold(1,j) + sum(v0tau(:,:,i)*v)
        Cvvold(2,j) =  Cvvold(2,j) + sum(v0s(:,:,i)*vs)
        Cvvold(3,j) =  Cvvold(3,j) + sum(p0(:,:,i)*v)
        j = j - 1
    end do

    if (naccumulated == window_width) then
        if (Cvvold(1,1) /= 0.0d0 ) then
            Cvvold = Cvvold / (Natom * window_width)
            call dump_cvv(Cvvold)
        end if
        Cvvold = Cvvcur
        Cvvcur = 0.0d0
        naccumulated = 0
    end if
end subroutine

subroutine mean_sq_disp(dxsq)
    real*8, intent(out) :: dxsq

    dxsq = sum((r - r0equil - rshift)**2)/Natom
end subroutine

subroutine dump_cvv(Cvv)
    use utils
    implicit none
    real*8, intent(in) :: Cvv (:,:)
    integer :: i
    character(4) :: cdump

    ncvvout = ncvvout + 1
    call int2strz(ncvvout, 4, cdump)
    open(cvvout, file=trim(stem)//'_Cvv_'//cdump//'.dat')
    write(cvvout,'(4F18.7)') (dt*(i-1)*t0fs, Cvv(1:3,i), i=1,size(Cvv,2))
    close(cvvout)
end subroutine
end module spine
