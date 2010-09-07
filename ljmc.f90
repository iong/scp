program pljmc
    use ljmc
    use mc
    use vgw
    use xyz
    implicit none
    include 'mpif.h'
    real*8 :: rcmin
    character(LEN=256) :: arg, inputf, fname
    integer :: i, j, n, NMC,mcblen,ierr
    character(4) :: csubme
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

    !call populate_cube(bl, rcmin, r(:,:))

    call int2str(me, csubme)
    fname = "r0."//csubme//".xyz"
    call load_xyz(r, fname)
    call vgw0(r(:,:), U0, beta(me+1), 0.0d0, y0)

    call mc_burnin(10000)
    !write (*,*) 'xstep', xstep
    
    !call int2str(me, csubme)
    !fname = "r0."//csubme//".xyz"
    !call dump_xyz(r, fname, "(pH2)_50")

    Z = 0.0d0
    do n=1,NMC,mcblen
        call mc_block(mcblen, 1000)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
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
contains

subroutine int2str(i, csubme)
    integer, intent(in) :: i
    character(4) :: csubme
    csubme =   achar(i/1000+48) &
           // achar(mod(i/100,10)+48) &
           // achar(mod(i/10,10)+48) &
           // achar(mod(i,10)+48)
end subroutine
end program pljmc
! vim:et
