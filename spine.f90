module spine
    implicit none
    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer, parameter :: eout = 30, cvvout = 31
    integer :: Natom, Nbath, ntracks, tracksep, seglen
    real*8 :: rcmin, tstart, tstop, dt, bl, bl2,rho, kT
    real*8, dimension(:), allocatable :: y, Qbath, xi, vxi
    real*8, dimension(:,:), allocatable :: r0,  r, p, v, rshift
    real*8, dimension(:,:,:), allocatable :: Qnk, Meff, invMeff, v0tau, v0s, p0, vkubo, track, r0corr, r0shift, r0k
    integer, allocatable :: trackstart(:)
    real*8 :: lastepot
    character(256) :: stem

    interface
        subroutine nose_hoover_chain(p, Ekin, kT, xi, vxi, Q, dt, ne)
            real*8, intent(inout) :: p(:,:), Ekin, xi(:), vxi(:)
            real*8, intent(in) :: kT, Q(:), dt
            integer, intent(in) :: ne
        end subroutine nose_hoover_chain
        subroutine kubo(q0, v, beta, nsteps, xk, vk)
            real*8, intent(in) :: q0(:,:), v(:,:), beta
            integer, intent(in) :: nsteps
            real*8, intent(out) :: xk(:,:), vk(:,:)
        end subroutine kubo
        subroutine correlations(ndt)
            integer, intent(in) :: ndt
        end subroutine correlations
        subroutine dump_track(track, trackno)
            real*8, intent(in) :: track (:,:)
            integer, intent(in) :: trackno
        end subroutine dump_track
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
    real*8 :: Ueff, fcm(3)
    integer :: i

    p = p + 0.5*dt*f
    do i=1,Natom
        r(:,i) = r(:,i) + dt*matmul(invMeff(:,:,i), p(:,i))
    end do
    call vgw1(r, Ueff, 1.0/kT, 0.0d0, y, Meff, invMeff)
    f = 2.0*kT*unpack_f(y)
    fcm = sum(f,2)/Natom
    do i=1,Natom
        f(:,i) = f(:,i) - fcm
    end do
    p = p + 0.5*dt*f

    epot = -2.0*kT*y(1)
end subroutine


real*8 function kinetic_energy() result(ekin)
    implicit none
    integer :: j

    ekin = 0.0d0
    do j=1,Natom
        ekin = ekin + dot_product(p(:,j), matmul(invMeff(:,:,j), p(:,j)))
    end do
    ekin = 0.5*ekin
end function
end module spine
