module vgw
        real*8, dimension(10) :: LJA, LJC
        integer :: NGAUSS
        integer, parameter :: MAXN = 5000
        logical QRC(MAXN*MAXN)
        real*8 :: MASS, BL, RC, ATOL, TAUMIN
        integer :: N_atom
end module vgw

subroutine vgwinit(imass, ng, c, a, boxlen, rc_, taumin_, atol_)
        use vgw
        implicit none
        integer, intent(in) :: ng
        real*8, intent(in) :: imass, boxlen, c(ng), a(ng), rc_, taumin_, atol_

	MASS=1.0D0/(0.020614788876D0*IMASS)
        NGAUSS = ng
        LJC(1:ng) = c(1:ng)
        LJA(1:ng) = a(1:ng)
	TAUMIN=taumin_
        BL = boxlen
	RC = rc_
	ATOL=atol_
end subroutine