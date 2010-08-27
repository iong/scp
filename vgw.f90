module vgw
        real*8, dimension(10) :: LJA, LJC
        integer :: NGAUSS
        logical, allocatable :: QRC(:)
        real*8 :: MASS, BL, RC, ATOL, RTOL, TAUMIN
        integer :: N_atom
        real*8 ::
contains
include 'potential_energy.f90'
end module vgw

subroutine vgwinit(natom, real_mass, ng, c, a, boxlen, rc_, taumin_, atol_, rtol_)
        use vgw
        implicit none
        integer, intent(in) :: natom, ng
        real*8, intent(in) :: real_mass, boxlen, c(ng), a(ng), rc_, taumin_, atol_

        N_atom = natom
        MASS=1.0D0/real_mass
        NGAUSS = ng
        LJC(1:ng) = c(1:ng)
        LJA(1:ng) = a(1:ng)
        TAUMIN=taumin_
        BL = boxlen
        RC = rc_
        ATOL=atol_
        RTOL=rtol_

        allocate(QRC(natom*natom))
end subroutine
