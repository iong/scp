module vgwspfm
    use integration
    use utils
    use sparse
    implicit none
    private
    public :: vgwspfminit,vgw0spfm,vgw0spfmgs,vgwspfmcleanup
    
    integer, parameter :: info = 1
    integer :: debug = 0
    integer :: Natom, NGAUSS, NEQ
    real*8 :: BL, Vcutoff, Gcutoff, LAMBDA, sigma0, epsilon0
    real*8, dimension(10) :: LJA, LJC

    type(csr), target :: G, UXY, GU
    type csr_list
        type(csr), pointer :: p
    end type csr_list

    real*8 :: gama, gamap, U
    real*8, allocatable :: Gbdiag(:,:), Gband(:,:), UXYbdiag(:,:), UXYband(:,:), UX(:)
    integer, allocatable :: Gband_ja(:)


    double precision, allocatable, target :: Y(:)
    class(integrator), pointer :: prop

    real*8 :: invmass, mass, dt0, dtmax, dtmin, vgw_atol(3)

contains

subroutine vgwspfminit(species, M, Vrange, Grange)
    implicit none
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, Vrange, Grange

    include 'species.f90'
end subroutine

subroutine vgwspfmcleanup()
     call G%cleanup()
     call UXY%cleanup()
     call GU%cleanup()
     call prop % cleanup()
     deallocate(Gbdiag, UXYbdiag, UX, y)
 end subroutine vgwspfmcleanup

include 'vgw0spfm.f90'

include 'rhssspfm.f90'
include 'rhssfm1.f90'

end module vgwspfm
