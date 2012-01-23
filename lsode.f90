module lsode_mod

    type :: lsode
        double precision :: dt0, dtmin, dtmax
        double precision, allocatable :: RWORK(:), YP(:), ATOL(:)
        integer, allocatable :: IWORK(:)

        integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
        double precision :: RTOL
    contains
        procedure :: init=> lsode_init
        procedure :: propagate => lsode_prop
        procedure :: get_ncalls => lsode_get_ncalls 
        procedure :: get_nsteps => lsode_get_nsteps 
        procedure :: destroy => lsode_destroy 
        procedure :: set_dt => lsode_set_dt
        procedure :: set_atol => lsode_set_atol
    end type lsode
contains

    subroutine lsode_init(this, NEQ)
        implicit none
        class(lsode) :: this
        integer :: NEQ

        this%LRW = 20 + 16*NEQ
        this%LIW = 30

        if (.not. allocated(this%ATOL)) then
            allocate(this%ATOL(NEQ), this%RWORK( this%LRW ), &
                this%IWORK( this%LIW ))
        end if

        this%ITOL=2
        this%RTOL=1d-5
        this%ITASK=1
        this%ISTATE=1
        this%IOPT = 1
        this%MF=10
        this%IWORK=0

        this%IWORK(6) = 50000 !MXSTEP

        this%RWORK(5) = this%dt0
        this%RWORK(6) = this%dtmax
        this%RWORK(7) = this%dtmin
    end subroutine


    subroutine lsode_set_atol(this, atol)
        implicit none
        class(lsode) :: this
        double precision, intent(in) :: atol(:)

        this%atol = atol
    end subroutine


    subroutine lsode_prop(this, F, y, T, TSTOP)
        implicit none
        class(lsode) :: this
        interface
            subroutine F(NEQ, T, Y, YP)
                integer, intent(in) :: NEQ
                double precision, intent(in) ::  T
                double precision, intent(in), target ::  Y(:)
                double precision, intent(out), target :: YP(:)
            end subroutine
        end interface
        double precision, intent(inout) :: y(:), T
        double precision, intent(in) :: TSTOP

        CALL DLSODE(F, size(y), Y, T, TSTOP, this%ITOL, this%RTOL, this%ATOL, &
            this%ITASK, this%ISTATE, this%IOPT,&
            this%RWORK, this%LRW, this%IWORK, this%LIW, JAC, this%MF)
    end subroutine

    subroutine lsode_destroy(this)
        implicit none
        class(lsode) :: this
        if (.not. allocated(this%ATOL) ) return

        deallocate(this%ATOL, this%RWORK, this%IWORK)
    end subroutine

    subroutine lsode_set_dt(this, dt0, dtmin, dtmax)
        implicit none
        class(lsode) :: this
        double precision , intent(in) :: dt0, dtmin, dtmax

        this%dt0 = dt0
        this%dtmin = dtmin
        this%dtmax = dtmax
    end subroutine

    function lsode_get_ncalls(this)
        implicit none
        class(lsode) :: this
        integer :: lsode_get_ncalls

        lsode_get_ncalls = this%IWORK(12)
    end function


    function lsode_get_nsteps(this)
        implicit none
        class(lsode) ::this
        integer :: lsode_get_nsteps

        lsode_get_nsteps = this%IWORK(11)
    end function


    subroutine JAC()
    end subroutine

end module lsode_mod
