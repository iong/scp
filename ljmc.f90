program pljmc
    use ljmc
    use mc
    use vgw
    implicit none
    include 'mpif.h'
    real*8 :: rcmin
    character(LEN=256) :: arg, inputf
    integer :: i, j, n, NMC,mcblen,ierr
    namelist /ljmccfg/Natom,imass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,NMC,mcblen,ptinterval,Tmin,Tmax,outfile,rho,rcmin

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    if (command_argument_count() == 0) then
        inputf='pH2.in'
    else if (command_argument_count() == 1) then
        call get_command_argument(1, arg)
        inputf=trim(arg)
    else
        write (*,*) 'too many arguments'
        stop
    endif


    call load_defaults()

    LJA=0.0d0
    LJC=0.0d0
    open(20, file=inputf)
    read(20, NML=ljmccfg)
    close(20)
    
    call setup_ljmc()
    call seed_rng()

    bl = (Natom/rho)**(1.0/3.0)
    bl2=bl/2

    if (rcmin**3 * rho >= 1.0) then
        write (*,"(A,F7.3)") 'rcmin is too large. Ensure that rcmin <', &
            rho**(-1.0/3.0)
        stop
    endif

    IMASS = IMASS*0.020614788876D0
    call vgwinit(natom, bl)
    call populate_cube(bl, rcmin, r(:,:))
    U0=0
    call vgw0(r(:,:), U0, beta(me+1), 0.0d0, y0)

    call mc_burnin(10)
    write (*,*) 'xstep', xstep


    Z = 0.0d0
    do n=1,NMC,mcblen
        call mc_block(mcblen, 1000)
        write (*,*) 'xstep', xstep
    
        call MPI_Allgather(Z(me+1), 1, MPI_REAL8, Z, 1, MPI_REAL8, MPI_COMM_WORLD, ierr)
        call heat_capacity(nprocs, Z, kT, Cv)
        if (me==0) then
            call dump_Tmin(n)
            open(30, file='Z.dat')
            write (30,'(4(ES16.8," "))') (kT(i), Cv(i), Z(i), beta(i),i=1,nprocs)
            close(30)
        endif
    enddo

    call dump_Tmin(NMC)
    call MPI_Finalize(ierr)
end program pljmc
! vim:et
