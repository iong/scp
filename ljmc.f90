program ljmc
    use ljmc
    real*8 :: taumin = 1.0e-6, rtol=1e-4, atol=1.0e-4, rc=10.0
    real*8, parameter, dimension(4) :: &
            LJA = (/1.038252215127D0, 0.5974039109464D0, &
                    0.196476572277834D0, 0.06668611771781D0/), &
            LJC = (/96609.488289873d0, 14584.62075507514d0, &
                    -365.460614956589d0, -19.5534697800036d0/)
    character, dimension(256) :: arg, inputf
    namelist /cfglist/Natom,imass,bl,rc,rtol,atol,taumin,NMC,nstreams,Tmin,Tmax,output,rhosolid,ncells

    if (command_argument_count() /= 1) then
        write (*,*) 'Need input file'
        stop
    endif

    call get_command_argument(1, arg)
    inputf=trim(arg)

    open(20, file=inputf)
    read(20, NML=cfglist)
    close(20)

    setup_ljmc()
    kT(1:nstreams) = (/(Tmin + (Tmax-Tmin)/(nstreams-1)*(i-1), i=1,nstreams)/)
    beta(1:nstreams) = 1.0/kT(1:nstreams)

    open(20, file='Ne147.dat')
    read(20, *) (q0(1,i), q0(2,i), q0(3,i),i=1,natom)
    close(20)

    call  vgwinit(natom, IMASS*0.020614788876D0, 4, LJC, LJA, bl, rc, taumin, atol, rtol)

    do n=1,NMC
        do i=1,nstreams
            call mc_1by1(mcburn, i)
            Z(i) = Z(i) + exp(-beta(i) * U0(i))
        enddo

        call swap_streams()

        if (mod(n,100) == 0) then
            call dump_Umin(i)
            call heat_capacity(Z, kT, Cv)
            open(30, file='Z.dat')
            write (30,'(4F16.8)') (kT(i), Cv(i), Z(i), beta(i),i=1,nstreams)
            close(30)
        endif
    enddo

    call dump_Umin(NMC)
end program ljmc
! vim:et
