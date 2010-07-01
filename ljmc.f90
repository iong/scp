program ljmc
        integer, parameter :: natom = 13
        real*8 :: q0(3,natom), qtau(3,natom), fx(3,natom), y(1+21*natom), invmeff(3,3,natom), sqrtmeff(3,3,natom)
        real*8, parameter :: imass = 2.0, taui = 1.0e-6, bl=40.0, atol=1.0e-4, rc=10.0
        real*8 :: lnp, w, enrg, taumax
        real*8, parameter, dimension(4) :: &
                LJA = (/1.03825, 0.597404, 0.1964765, 0.0666861 /), &
                LJC = (/96609.488289873, 14584.6207550751, &
                        -365.460614956589, -19.5534697800036/)


        open(20, file='Ne13.dat')
        read(20, *) (q0(1,i), q0(2,i), q0(3,i),i=1,natom)
        close(20)

        taumax = 1.0

        call  vgwinit(bl, 4, LJC, LJA)
        qtau = q0
        call vgw0(natom,IMASS,qtau,W, TAUMAX,TAUI,ATOL,RC,Y)
        write (*,*) 'Ueff =', W
        qtau = q0
        call vgwquenchspb(natom,IMASS,qtau,FX,LNP,W, ENRG,TAUMAX,TAUI,BL,ATOL,RC,Y, InvMeff,SqrtMeff)
        write (*,*) 'Ueff =', W
        qtau = q0
        write (*,*) 'Ueff =', W
end program ljmc
! vim:et
