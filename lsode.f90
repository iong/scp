module lsode_mod
    use integrator_mod
    private

    type, public, extends(integrator) ::  lsode
        integer :: LIW, MF
        integer :: IOPT,ISTATE,LRW,ITOL,IERR
        integer, allocatable :: IWORK(:)
        double precision, allocatable :: RWORK(:)
    contains
        procedure :: init
        procedure :: advance
        procedure :: get_xp
        procedure :: forget_history
        procedure :: set_dt
        procedure :: set_dtmin
        procedure :: set_dtmax
        !procedure :: step => lsode_step
        procedure :: cleanup
    end type lsode
contains
    subroutine  init(self, NEQ)
        implicit none
        class(lsode) :: self
        integer, intent(IN) :: NEQ

        call self % integrator % init(NEQ)

        self % LRW = 20 + 16*NEQ

        self % LIW = 30
        self % MF = 10
        self % ITOL=2
        self % IOPT=1
        self % ISTATE=1

        allocate(self % RWORK( self % LRW), self % IWORK(self % LIW))

        self % IWORK=0
        self % IWORK(6) = 50000 ! MXSTEP

        self % RWORK(5:10)=0.0D0
    end subroutine init


    subroutine set_dt(self, dt)
        class(lsode) :: self
        double precision, intent(in) :: dt

        self % dt = dt
        self%RWORK(5) = self % dt
    end subroutine set_dt


    subroutine set_dtmin(self, dtmin)
        class(lsode) :: self
        double precision, intent(in) :: dtmin

        self % dtmin = dtmin
        self%RWORK(7) = self % dtmin
    end subroutine set_dtmin



    subroutine set_dtmax(self, dtmax)
        class(lsode) :: self
        double precision, intent(in) :: dtmax

        self % dtmax = dtmax
        self%RWORK(6) = self % dtmax
    end subroutine set_dtmax


    subroutine forget_history(self)
        class(lsode) :: self

        self % ISTATE = 1
    end subroutine


    subroutine advance(self, F, x, tstop)
        implicit none
        class(lsode) :: self
        DOUBLE PRECISION, intent(inout) :: x(:)
        double precision, intent(in) :: tstop
        procedure(RHS_X) :: F

        CALL DLSODE(f77_rhs, size(x), x, self % t, tstop, self % ITOL, &
                self % RTOL, &
                self % ATOL, self % TASK, self % ISTATE, self % IOPT, &
                self % RWORK, self % LRW, self % IWORK, self % LIW, JAC, &
                self % MF)

        self % dt = self % RWORK(12)
        self % nsteps = self % IWORK(11)
        self % ncalls = self % IWORK(12)

    contains
        subroutine f77_rhs(neq, t, y, yp)
            integer, intent(in) :: NEQ
            double precision, intent(in) :: T
            double precision, intent(in), target :: Y(NEQ)
            double precision, intent(out), target :: YP(NEQ)

            call F(neq, t, y, yp)
        end subroutine f77_rhs
    end subroutine advance


    subroutine get_xp(self, xp)
        implicit none
        class(lsode) :: self
        DOUBLE PRECISION, intent(out) :: xp(:)
        integer :: iflag

        call DINTDY(self%t, 1, self%RWORK(21), self%neq, xp, iflag)
    end subroutine


    subroutine cleanup(self)
        class(lsode) :: self

        deallocate(self % RWORK, self % IWORK)

        call self % integrator % cleanup()
    end subroutine cleanup


    subroutine  JAC()
    end subroutine JAC


end module
