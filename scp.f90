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

    double precision :: dkT=0.03d0, kTstop=0.3, kT, Gcutoff=2.75, Vcutoff=2.75

    double precision, allocatable :: r0(:,:)
    double precision :: Uref, E0

    class(vgw), pointer :: p

    nargs = command_argument_count()
    narg = 0

    do
        narg = narg + 1
        if (narg > nargs) exit
        call get_command_argument(narg, arg)
        if (trim(arg) == '-sparse') then
            allocate(ffvgw::p)
            cycle
        end if

        if (trim(arg) == '-rc') then
            call get_command_argument(narg+1, arg)
            read(arg, *) Gcutoff
            call get_command_argument(narg+2, arg)
            read(arg, *) Vcutoff

            narg = narg + 2
            cycle
        end if

        call load_xyz(arg, r0)
       
        call get_command_argument(narg+1, arg)
        read(arg, *) dkT

        call get_command_argument(narg+2, arg)
        read(arg, *) kTstop
        narg = narg + 2
    end do

    Natom = size(r0, 2)

    if (.NOT. associated(p)) allocate(vgwfm::p)

    call p%init(Natom, 'LJ')
    ! set the cutoffs if the sparse approach is used

    select type (p)
    type is(ffvgw)
        call p%set_range(Gcutoff, Vcutoff)
    end select

    kT = 0d0
    E0 = p%Utot0(r0)/Natom
    print '(F6.3,2F12.7,2E12.5)', kT, E0, 0d0, 0d0

    do i=1,nint(kTstop/dkT)
        kT = dkT*i
        E0 = p%F(reshape(r0, (/3*Natom/)) , kT, kT>dkT)
        print '(F6.3,2F12.7,2E12.5)', kT, E0, p % qconv, p % gconv
    end do

    write(waste, '("rt",I,".dat")') Natom
    open(38,file=trim(waste))
    write(38, '(I10,2F16.7,I10)') Natom, p % iter_time / p % niterations, &
          p % logdet_time * kT / dkT,  p % niterations
    call p%cleanup()
    deallocate(p)

end program scp
