program pljmc
    use ljmc
    implicit none
    real*8 :: taumin = 1.0e-6, rtol=1e-4, atol=1.0e-4, rc=10.0
    real*8, parameter, dimension(4) :: &
            LJA = (/1.038252215127D0, 0.5974039109464D0, &
                    0.196476572277834D0, 0.06668611771781D0/), &
            LJC = (/96609.488289873d0, 14584.62075507514d0, &
                    -365.460614956589d0, -19.5534697800036d0/)
    real*8 :: rcmin, rcminsq, dsq
    character(LEN=256) :: arg, inputf
    logical :: too_close
    integer :: i, j, n, NMC, mcburn
    namelist /ljmccfg/Natom,imass,bl,rc,rtol,atol,taumin,NMC,mcburn,nstreams,Tmin,Tmax,outfile,rho,rcmin,ncells

    if (command_argument_count() /= 1) then
        write (*,*) 'Need input file'
        stop
    endif

    call get_command_argument(1, arg)
    inputf=trim(arg)

    open(20, file=inputf)
    read(20, NML=ljmccfg)
    close(20)

    call setup_ljmc()

!    open(20, file='Ne147.dat')
!    read(20, *) (q0(1,i), q0(2,i), q0(3,i),i=1,natom)
!    close(20)

    kT(1:nstreams) = (/(Tmin + (Tmax-Tmin)/(nstreams-1)*(i-1), i=1,nstreams)/)
    beta(1:nstreams) = 1.0/kT(1:nstreams)
    bl = (Natom/rho)**(1.0/3.0)
    bl2=bl/2

    rcminsq = rcmin*rcmin
    if (rcminsq*rcmin*rho >= 1.0) then
        write (*,"(A,F7.3)") 'rcmin is too large. Ensure that rcmin <', &
            rho**(-1.0/3.0)
        stop
    endif

    do i=1,Natom
        do
            too_close = .FALSE.
            call random_number(q0(:,i))
            q0(:,i) = q0(:,i) * bl
            do j=1,i-1
                dsq = sum((q0(:,j)-q0(:,i))**2)
                if (dsq < rcminsq) then
                    too_close = .TRUE.
                    exit
                endif
            enddo
            if (.not. too_close) exit
        enddo
    enddo
    r = spread(q0, 3, nstreams)

    call  vgwinit(natom, IMASS*0.020614788876D0, 4, LJC, LJA, bl, rc, taumin, atol, rtol)

    do n=1,NMC
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(r,rnew,rmin,U0,Umin,rmove,Z,kT,beta,naccepted,ntrials)
!$OMP DO SCHEDULE(GUIDED) 
        do i=1,nstreams
            call mc_1by1(mcburn, i)
            Z(i) = Z(i) + exp(-beta(i) * U0(i))
        enddo
!$OMP END DO
!$OMP END PARALLEL

        call pt_swap()

        if (mod(n,2) == 0) then
            call dump_Umin(i)
            call heat_capacity(Natom, Z, kT, Cv)
            open(30, file='Z.dat')
            write (30,'(4F16.8)') (kT(i), Cv(i), Z(i), beta(i),i=1,nstreams)
            close(30)
        endif
    enddo

    call dump_Umin(NMC)
end program pljmc
! vim:et
