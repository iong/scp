module vgw
        real*8, dimension(10) :: LJA, LJC
        integer :: NGAUSS
        integer, parameter :: MAXN = 5000
        logical QRC(MAXN*MAXN)
        real*8 :: MASS, BL
        integer :: N_atom
end module vgw

subroutine vgwinit(boxlen, ng, c, a)
        use vgw
        implicit none
        integer, intent(in) :: ng
        real*8, intent(in) :: boxlen, c(ng), a(ng)

        LJA(1:ng) = a(1:ng)
        LJC(1:ng) = c(1:ng)
        NGAUSS = ng
        BL = boxlen
end subroutine
