program clustergs
!    use vgw
    use ffvgw_mod
    use vgwfm_mod
    use utils, only: M_PI
    implicit none
!    include 'mkl_service.fi'

    character(256) :: arg, coords, waste

    integer :: i, Natom, status

    double precision :: kT=0.02d0, kTstop=0.3

    double precision, allocatable :: r0(:,:)
    double precision :: Uref, E0

    class(ffvgw), pointer :: fm

    allocate(fm)

    call get_command_argument(1, coords)

    open(33, file=trim(coords))
    read (33, *) Natom

    allocate(r0(3,Natom))

    Uref = 0d0
    read(33,IOSTAT=status,FMT=*) Uref
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
    if (command_argument_count() > 2) then
        call get_command_argument(3, arg)
        read(arg, *) kTstop
    end if

    call fm%init(Natom, 'LJ')
    call fm%set_range(1.8d0, 2.75d0)
    E0 = fm%F(reshape(r0, (/3*Natom/)) , kT)

    !fm % Omega % x => fm%y(3*Natom+1 : 3*Natom + fm%Omega%nnz)
    !call fm % Omega % write('Omega.csr')
    !nullify(fm % Omega % x)

    do
        print '(F6.3,2F12.7,2E12.5)', kT, E0, fm % qconv, fm % gconv
        kT = kT + 0.01
        if (kT >  kTstop) exit
        E0 = fm%F(reshape(r0, (/3*Natom/)) , kT, .TRUE.)
    end do
    call fm%cleanup()

end program clustergs
