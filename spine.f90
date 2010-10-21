module spine
    implicit none
    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer :: Natom, Nbath
    real*8 :: rcmin, tstart, tstop, dtout, dt, bl, bl2,rho, kT
    real*8, dimension(:), allocatable :: y, Qbath, xi, vxi
    real*8, dimension(:,:), allocatable :: r0, p0, vtau0, v0, r(:,:), p(:,:), rshift(:,:), r0equil(:,:)
    real*8, dimension(:,:,:), allocatable :: Meff0, invMeff0, Qnk0, Qnk, Meff, invMeff
    real*8 :: lastepot

    interface
        subroutine nose_hoover_chain(p, Ekin, kT, xi, vxi, Q, dt, ne)
            real*8, intent(inout) :: p(:,:), Ekin, xi(:), vxi(:)
            real*8, intent(in) :: kT, Q(:), dt
            integer, intent(in) :: ne
        end subroutine nose_hoover_chain
    end interface
contains


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
    real*8 :: Ueff, v3(3), p3(3), Ekin
    integer :: i, Nq, Nqp


    Natom = size(r, 2)
    p = p + 0.5*dt*f
    do i=1,Natom
        r(:,i) = r(:,i) + dt*matmul(invMeff(:,:,i), p(:,i))
    end do
    call vgw1(r, Ueff, 1.0/kT, 0.0d0, y, Meff, invMeff)
    f = reshape( 2.0*kT*y(2+18*Natom:1+21*Natom), (/ 3, Natom /))
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

subroutine velocity_autocorrelation(Cvv)
    implicit none
    real*8, intent(out) :: Cvv
    real*8 :: v3(3)
    integer :: j

    Cvv = 0.0d0
    do j=1,Natom
        v3 = matmul(invMeff(:,:,j), p(:,j))
        Cvv = Cvv + dot_product(vtau0(:,j), v3)
    end do
    Cvv = Cvv / Natom
end subroutine

subroutine mean_sq_disp(dxsq)
    real*8, intent(out) :: dxsq

    dxsq = sum(sum((r - r0equil - rshift)**2, 1))/Natom
end subroutine
end module spine
