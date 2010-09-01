module ljmc
    implicit none
    integer :: Natom
    real*8 :: M_PI = 3.14159265358979323846264338327950288d0
    real*8 :: bl=1d10, bl2,rho
    integer::ncells
    real*8,allocatable :: q0(:,:), y0(:)
    character(len=256) :: outfile

    interface
    real*8 function gaussran(sigma, x0) result(y)
        real*8, intent(in) :: sigma, x0
    end function
    end interface

    interface
    subroutine heat_capacity(N, Z, kT, Cv)
        integer, intent(in) :: N
        real*8, intent(in), dimension(N) :: Z, kT
        real*8, intent(out), dimension(N) :: Cv
    end subroutine
    end interface

    interface
    subroutine populate_cube(bl, rcmin, r)
        real*8, intent(in) :: bl, rcmin
        real*8, intent(out) :: r(:,:)
    end subroutine
    end interface

    interface
    subroutine load_defaults()
    end subroutine
    end interface
end module ljmc
