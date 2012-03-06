program scp
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

    double precision :: kT=0.03d0, kTstop=0.3

    double precision, allocatable :: r0(:,:)
    double precision :: Uref, E0

    class(vgw), pointer :: p

    nargs = command_argument_count()
    narg = 1

    call get_command_argument(narg, arg)
    if (trim(arg) == '-sparse') then
        allocate(ffvgw::p)
        narg = narg + 1
    else
        allocate(vgwfm::p)
    end if
        
    call get_command_argument(narg, arg)
    call load_xyz(arg, r0)
    narg = narg + 1
    Natom = size(r0, 2)

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

    call p%init(Natom, 'LJ')
    ! set the cutoffs if the sparse approach is used
    select type (p)
    type is(ffvgw)
        call p%set_range(1.8d0, 2.75d0)
    end select

    ! kT = 0
    E0 = p%Utot0(r0)/Natom
    print '(F6.3,2F12.7,2E12.5)', 0d0, E0, 0d0, 0d0

    E0 = p%F(reshape(r0, (/3*Natom/)) , kT)

    !ff % Omega % x => ff%y(3*Natom+1 : 3*Natom + ff%Omega%nnz)
    !call ff % Omega % write('Omega.csr')
    !nullify(ff % Omega % x)

    do
        print '(F6.3,2F12.7,2E12.5)', kT, E0, p % qconv, p % gconv
        kT = kT + 0.03
        if (kT >  kTstop) exit
        E0 = p%F(reshape(r0, (/3*Natom/)) , kT, .TRUE.)
    end do
    call p%cleanup()

end program scp
