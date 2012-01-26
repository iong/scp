
!!!!!!!
!!!!!!! solve damped harmonic oscillator in q, G => this
!!!!!!! will move q and G towards the minimum of the potential energy F (free energy)


program clustergs
!    use vgw
    use vgwfm_mod
    use utils, only: M_PI
    implicit none
!    include 'mkl_service.fi'

    character(256) :: arg, coords, waste

    integer :: i, Natom, status

    double precision :: endtime(3), begtime(3), mass, taustop
    double precision ::  deBoer=0.1, kT=0.001d0, RC=100d0, BL=100d0, sigma0 = 2.749, epsilon0 = 35.6d0

    double precision, allocatable :: r0(:,:)
    double precision :: Uref, E0

    class(vgwfm), pointer :: fm

    allocate(fm)

    call get_command_argument(1, coords)

    open(33, file=trim(coords))
    read (33, *) Natom

    allocate(r0(3,Natom))

    Uref = 0d0
    read(33,IOSTAT=status,FMT=*) deBoer, Uref
    if (status > 0) then
        Uref = 0d0
    else if (status<0) then
        write (*,*) 'read(33,*) failed with status', status
        stop
    end if

    read(33, *) (waste, r0(:,i), i=1,Natom)
    close(33)

    ! deBoer can be supplied as a second argument on the command line
    if (command_argument_count() > 1) then
        call get_command_argument(2, arg)
        read(arg, *) deBoer
    end if

    r0 = r0*1.05d0 / sigma0

    taustop = 1d0/kT
    mass = 1/deBoer**2

    call fm%init(Natom, 'LJ')
    call fm%set_mass((/ ( mass, i=1,Natom ) /) )
    E0 = fm%F(r0, kT, taustop)
    call fm%cleanup()

    !U(2:4) = U(2:4) * epsilon0

    write(*,'(I10,2F12.7)') Natom, E0/Natom, sqrt(sum(r0**2)/(3*Natom))

end program clustergs
