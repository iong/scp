module vgw
    real*8, dimension(10), private :: LJA, LJC
    integer, private :: NGAUSS
    real*8, private :: BL
    logical, allocatable :: QRC(:)
    real*8 :: MASS, RC, ATOL, RTOL, TAUMIN, IMASS
    integer :: N_atom
    
contains

subroutine vgwinit(natom, ng, c, a, boxlen)
    implicit none
    integer, intent(in) :: natom, ng
    real*8, intent(in) :: boxlen, c(ng), a(ng)
    N_atom = natom
    NGAUSS = ng
    LJC(1:ng) = c(1:ng)
    LJA(1:ng) = a(1:ng)
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

