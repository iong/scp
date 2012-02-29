program clustergs
!    use vgw
    use ffvgw_mod
    use vgwfm_mod
    use utils, only: M_PI
    implicit none
!    include 'mkl_service.fi'

    character(256) :: arg, coords, waste

    integer :: i, Natom, status

    double precision :: kT=0.001d0

    double precision, allocatable :: r0(:,:)
    double precision :: Uref, E0

    class(ffvgw), pointer :: fm

    allocate(fm)

    call get_command_argument(1, coords)

    open(33, file=trim(coords))
    read (33, *) Natom

    allocate(r0(3,Natom))

    Uref = 0d0
    read(33,IOSTAT=status,FMT=*) kT, Uref
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
        read(arg, *) kT
    end if

    r0 = r0! / sigma0

    call fm%init(Natom, 'LJ')
    call fm%set_range(5d0, 5d0)
    E0 = fm%F(reshape(r0, (/3*Natom/)) , kT)
    call fm%cleanup()

    !U(2:4) = U(2:4) * epsilon0

    write(*,'(I10,2F12.7)') Natom, E0, fm % gconv

end program clustergs
