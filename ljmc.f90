program pljmc
    use ljmc
    use mc
    use vgw
    use xyz
    use utils
    implicit none
    include 'mpif.h'
    real*8 :: rcmin
    character(LEN=256) :: arg, inputf, fname
    integer :: i, j, n, NMC,mcblen,ierr,mcburn
    character(4) :: csubme
    namelist /ljmccfg/Natom,imass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,NMC,mcblen,mcburn,ptinterval,Tmin,Tmax,outfile,rho,rcmin

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

    call int2strz(me, 4, csubme)
    fname = "r0."//csubme//".xyz"
    call load_xyz(r, fname)
    call vgw0(r(:,:), U0, taugrid, 0.0d0, y0)

    call mc_burnin(mcburn)

    call mc_run_loop(NMC, mcblen, 1000)

    call MPI_Finalize(ierr)
end program pljmc
! vim:et
