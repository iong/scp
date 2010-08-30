module ljmc
    integer :: Natom
    real*8 :: M_PI = 3.14159265358979323846264338327950288d0
    real*8 :: bl=1e10, bl2, Tmin, Tmax,rho,imass
    integer::ncells, nstreams
    real*8,allocatable :: q0(:,:), y0(:)
    real*8, allocatable :: r(:,:,:), rnew(:,:,:), rmin(:,:,:), U0(:), Umin(:)
    real*8, allocatable :: rmove(:), Z(:), kT(:), beta(:), Cv(:)
    integer, allocatable :: naccepted(:), ntrials(:)
    character(len=256) :: outfile
contains
real*8 function gaussran(sigma, x0) result(y)
    implicit none
    real*8, intent(in) :: sigma, x0
    real*8 :: x(2), M_PI = 3.14159265358979323846264338327950288
    call random_number(x)
    do while (x(1) == 0)
        call random_number(x)
    enddo
    y = sqrt( -2.0 * log(x(1))) * cos(2*M_PI*x(2))

    y = y*sigma + x0
end function

subroutine setup_ljmc()
    allocate (q0(3,natom), y0(1+21*natom), &
            r(3,natom,nstreams), rnew(3,natom,nstreams), &
            rmin(3,natom,nstreams), U0(nstreams), &
            Umin(nstreams), ntrials(nstreams), naccepted(nstreams),rmove(nstreams),&
            Z(nstreams), kT(nstreams), beta(nstreams), Cv(nstreams))
end subroutine
end module ljmc

