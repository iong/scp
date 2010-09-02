program pljmc
    use ljmc
    use mc
    use vgw
    implicit none
    include 'mpif.h'
    real*8, parameter, dimension(4) :: &
            LJA = (/1.038252215127D0, 0.5974039109464D0, &
                    0.196476572277834D0, 0.06668611771781D0/), &
            LJC = (/96609.488289873d0, 14584.62075507514d0, &
                    -365.460614956589d0, -19.5534697800036d0/)
    real*8, parameter, dimension(3) :: &
            pH2A = (/0.680688416015581, 0.172533385762719, 0.0509575639089336/), &
            pH2C = (/31319.4442342458, -282.50951658104, -9.27807987319189/)
            
    real*8 :: rcmin
    character(LEN=256) :: arg, inputf
    integer :: i, j, n, NMC,mcblen,ierr
    namelist /ljmccfg/Natom,imass,rc,rtol,atol,taumin,NMC,mcblen,ptinterval,Tmin,Tmax,bl,outfile,rho,rcmin

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    if (command_argument_count() /= 1) then
        write (*,*) 'Need input file'
        stop
    endif

    call get_command_argument(1, arg)
    inputf=trim(arg)

    call load_defaults()

    open(20, file=inputf)
    read(20, NML=ljmccfg)
    close(20)
    
    call setup_ljmc()
    call seed_rng()

!    open(20, file='Ne147.dat')
!    read(20, *) (q0(1,i), q0(2,i), q0(3,i),i=1,natom)
!    close(20)

    bl = (Natom/rho)**(1.0/3.0)
    bl2=bl/2

    if (rcmin**3 * rho >= 1.0) then
        write (*,"(A,F7.3)") 'rcmin is too large. Ensure that rcmin <', &
            rho**(-1.0/3.0)
        stop
    endif

    IMASS = IMASS*0.020614788876D0
    call vgwinit(natom, 3, pH2C, pH2A, bl)
    call populate_cube(bl, rcmin, r(:,:))
    call vgw0(r(:,:), U0, beta(me), 0.0d0, y0)

    call mc_burnin(10000)
    write (*,*) 'xstep', xstep


    do n=1,NMC,mcblen
        call mc_block(mcblen, 1000)
        write (*,*) 'xstep', xstep
    
        call MPI_Gather(Z(me), 1, MPI_REAL8, Z, nprocs, MPI_REAL8, 0, MPI_COMM_WORLD)
        if (me==0) then
            call dump_Tmin(n)
            call heat_capacity(nprocs, Z, kT, Cv)
            open(30, file='Z.dat')
            write (30,'(4(ES16.8," "))') (kT(i), Cv(i), Z(i), beta(i),i=1,nprocs)
            close(30)
        endif
    enddo

    call dump_Tmin(NMC)
    call MPI_Finalize(ierr)
end program pljmc
! vim:et
