program ljmc
        integer, parameter :: natom = 13
        real*8 :: q0(3,natom), fx(3,natom), y(1+21*natom), invmeff(3,3,natom), sqrtmeff(3,3,natom)
        real*8, parameter :: imass = 2.0, taui = 1.0e-6, bl=40.0, atol=1.0e-4, rc=10.0
        real*8 :: lnp, w, enrg, taumax

        open(20, file='Ne13.dat')
        read(20, *) (q0(1,i), q0(2,i), q0(3,i),i=1,natom)
        close(20)

        taumax = 1.0/20.0
        call vgwquenchspb(natom,IMASS,q0,FX,LNP,W, ENRG,TAUMAX,TAUI,BL,ATOL,RC,Y, InvMeff,SqrtMeff)
end program ljmc
