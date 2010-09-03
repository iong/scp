module vgw
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    real*8, private :: BL
    logical, allocatable :: QRC(:)
    real*8 :: MASS, RC, ATOL, RTOL, TAUMIN, IMASS
    integer :: N_atom
    
contains

subroutine vgwinit(natom, boxlen)
    implicit none
    integer, intent(in) :: natom
    real*8, intent(in) :: boxlen
    N_atom = natom
    BL = boxlen
    MASS=1.0D0/IMASS

    allocate(QRC(natom*natom))
end subroutine

include 'potential_energy.f90'
include 'vgw0.f90'
include 'vgw1.f90'
include 'rhss0.f90'
include 'rhss1.f90'

end module vgw

