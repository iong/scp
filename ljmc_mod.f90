module ljmc
    implicit none
    integer :: Natom
    integer :: me, nprocs
    real*8 :: bl=1d10, bl2,rho
    real*8,allocatable :: q0(:,:), y0(:)
    character(len=256) :: outfile

    interface
    subroutine heat_capacity(Z, kT, Cv)
        real*8, intent(in), dimension(:) :: Z, kT
        real*8, intent(out), dimension(:) :: Cv
    end subroutine

    subroutine heat_capacity2(Z, beta, Cv)
        real*8, intent(in), dimension(:) :: Z, beta
        real*8, intent(out), dimension(:) :: Cv
    end subroutine

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
