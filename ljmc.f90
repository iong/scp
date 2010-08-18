program ljmc
        integer, parameter :: natom = 147, ntests=10
        real*8 :: q0(3,natom), dq0(3,natom), q00(3,natom), y0(1+18*natom), fx(3,natom), y(1+21*natom), invmeff(3,3,natom), sqrtmeff(3,3,natom)
        real*8, parameter :: imass = 2.0, taui = 1.0e-6, bl=40.0, atol=1.0e-4, rc=10.0
        real*8 :: lnp, w, enrg, taumax, et(4)
        real*8, parameter, dimension(4) :: &
                LJA = (/1.03825, 0.597404, 0.1964765, 0.0666861 /), &
                LJC = (/96609.488289873, 14584.6207550751, &
                        -365.460614956589, -19.5534697800036/)
        real*8, dimension(ntests,2) :: Ueff
        integer :: seed_size, seed(10)
!POMP$ INST INIT 
        open(20, file='Ne147.dat')
        read(20, *) (q0(1,i), q0(2,i), q0(3,i),i=1,natom)
        close(20)

        taumax = 1.0

        call random_seed(SIZE=seed_size)
!        write(*,*) 'Seed size =', seed_size
        call random_seed(GET=seed(1:seed_size))
        call random_seed(PUT=seed(1:seed_size))

        call  vgwinit(natom, IMASS/48.5086, 4, LJC, LJA, bl, rc, taui, atol)
        q00 = q0
       call cpu_time(et(1))
!$OMP PARALLEL SHARED(Y, UPV_SHARED, UPM_SHARED) PRIVATE(tid,nthreads) FIRSTPRIVATE(NEQ,ISTATE,T)
        do i=1,ntests
                call random_number(dq0)
                q0 = q00 + dq0*0.1
                call vgw0(q0,Ueff(i, 1), TAUMAX,0.0, y0)
        enddo
!$OMP END PARALLEL
        call cpu_time(et(2))
!$POMP INST BEGIN(LJMCSTOP)
        call random_seed(PUT=seed(1:seed_size))
        call cpu_time(et(3))
        do i=1,ntests
                call random_number(dq0)
                q0 = q00 + dq0*0.1
                call vgwquenchspb(natom,IMASS,q0,FX,LNP,Ueff(i,2), ENRG,TAUMAX,TAUI,BL,ATOL,RC,Y, InvMeff,SqrtMeff)
        enddo
        call cpu_time(et(4))
        write (*,*) et(2)-et(1), et(4)-et(3)
        write(*,"(2F18.8)") ((Ueff(i, j),j=1,2),i=1,ntests)
!$POMP INST END(LJMCSTOP)
end program ljmc
! vim:et
