module spine
    implicit none
    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer, parameter :: eout = 30, cvvout = 31, track_width = 14
    integer :: Natom, ndt
    real*8 :: rcmin, tstart, tstop, dt, bl, bl2,rho, kT

    real*8, dimension(:,:), allocatable :: r0
    real*8, dimension(:), allocatable :: y
    real*8, dimension(:,:), allocatable :: r, p0, p, v, rshift, f, trackaccum, rkold
    real*8, dimension(:,:), allocatable :: v0tau, v0s, v0, vkubo, track, q0tau, r0shift, r0k, rprev, vprev, vsprev, r0s, rsprev,v0k
    real*8, dimension(:,:,:), allocatable :: Qnk, Meff, invMeff

    real*8 :: lastepot
    character(256) :: stem

    interface
        subroutine kubo(q0, v, beta, nsteps, xk, vk)
            real*8, intent(in) :: q0(:,:), v(:,:), beta
            integer, intent(in) :: nsteps
            real*8, intent(out) :: xk(:,:), vk(:,:)
        end subroutine kubo
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


subroutine verletstep(dt, Epot)
    use vgw
    implicit none
    real*8, intent(in) :: dt
    real*8, intent(out) :: Epot
    real*8 :: Ueff, fcm(3)
    integer :: i, k

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

    Epot = -2.0*kT*y(1)
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
