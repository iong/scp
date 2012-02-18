module ffvgw_mod
    use integration
    use utils
    use sparse
    use vgw_mod
    implicit none
    private

    type, public, extends(vgw) :: ffvgw
	double precision :: rfullmatsq
        double precision, allocatable :: UX(:)

        double precision :: Vcutoff, Gcutoff

        type(csr), target :: G, UXY, GU
        type csr_list
            type(csr), pointer :: p
        end type csr_list

        double precision, allocatable :: Gbdiag(:,:), Gband(:,:), UXYbdiag(:,:), UXYband(:,:), UX(:)
        integer, allocatable :: Gband_ja(:)
    contains
	procedure :: cleanup
        procedure :: converge
        procedure :: set_range
	procedure :: init_prop
	procedure :: logdet
end type ffvgw
    
contains

    subroutine set_range(self, Vcutoff, Gcutoff)
        implicit none
        class(ffvgw) :: self
        double precision :: intent(in) :: Vrange, Grange

        self % Vcutoff = Vcutoff
        self % Gcutoff = Gcutoff
    end subroutine set_range

    subroutine cleanup(self)
        class(ffvgw) :: self
        call self % G   % cleanup()
        call self % UXY % cleanup()
        call self % GU  % cleanup()

        deallocate(Gbdiag, UXYbdiag, UX, y)
    end subroutine cleanup

include 'vgw0spfm.f90'

include 'rhssspfm.f90'
include 'rhssfm1.f90'

end module vgwspfm
