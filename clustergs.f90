program clustergs
!    use vgw
    use vgw_mod
    use ffvgw_mod
    use vgwfm_mod
    use utils, only: M_PI
    use xyz, only: load_xyz
    implicit none
!    include 'mkl_service.fi'

    character(256) :: arg, coords, waste

    integer :: i, Natom, status, narg, nargs

    double precision :: kT=0.02d0, kTstop=0.3

    double precision, allocatable :: r0(:,:)
    double precision :: Uref, E0

    class(vgw), pointer :: scp

    nargs = command_argument_count()
    narg = 1

    call get_command_argument(narg, arg)
    if (arg == '-sparse') then
        allocate(ffvgw::scp)
        narg = narg + 1
    else
        allocate(vgwfm::scp)
    end if
        
    call get_command_argument(narg, arg)
    call load_xyz(arg, r0)
    narg = narg + 1

    if (narg <= nargs) then
        call get_command_argument(narg, arg)
        read(arg, *) kT
        narg = narg + 1
    end if

    if (narg <= nargs) then
        call get_command_argument(narg, arg)
        read(arg, *) kTstop
        narg = narg + 1
    end if

    call scp%init(Natom, 'LJ')
    ! set the cutoffs if the sparse approach is used
    select type (scp)
    type is(ffvgw)
        call scp%set_range(1.8d0, 2.75d0)
    end select
    E0 = scp%F(reshape(r0, (/3*Natom/)) , kT)

    !ff % Omega % x => ff%y(3*Natom+1 : 3*Natom + ff%Omega%nnz)
    !call ff % Omega % write('Omega.csr')
    !nullify(ff % Omega % x)

    do
        print '(F6.3,2F12.7,2E12.5)', kT, E0, scp % qconv, scp % gconv
        kT = kT + 0.01
        if (kT >  kTstop) exit
        E0 = scp%F(reshape(r0, (/3*Natom/)) , kT, .TRUE.)
    end do
    call scp%cleanup()

end program clustergs
