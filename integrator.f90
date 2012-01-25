module integrator_mod
    private

    type, public :: integrator
        double precision :: t, dt, dtmin, dtmax, rtol
        double precision, allocatable :: atol(:)

        integer :: nsteps, ncalls
    contains
        procedure :: init
        procedure :: set_dtmin
        procedure :: set_dtmax
        procedure :: advance
        procedure :: cleanup
    end type integrator
    
    abstract interface
        subroutine RHS_X(neq, t, y, yp)
            integer, intent(in) :: neq
            double precision, intent(in) :: t
            double precision, intent(in), target :: y(:)
            double precision, intent(out), target :: yp(:)
        end subroutine RHS_X
    end interface

    public :: RHS_X
contains

    subroutine init(self, NEQ, tstart, dt0)
        implicit none
        class(integrator) :: self
        integer, intent(IN) :: NEQ
        double precision, intent(in) :: tstart, dt0
                
        self % t = tstart
        self % dt = dt0
        self % dtmin = 0d0
        self % dtmax = 0d0

        self % nsteps = 0
        self % ncalls = 0

        allocate(self % atol(NEQ))
    end subroutine init


    subroutine set_dtmin(self, dtmin)
        class(integrator) :: self
        double precision, intent(in) :: dtmin

        self % dtmin = dtmin
    end subroutine set_dtmin


    subroutine set_dtmax(self, dtmax)
        class(integrator) :: self
        double precision, intent(in) :: dtmax

        self % dtmax = dtmax
    end subroutine set_dtmax


    subroutine advance(self, F, x, tstop)
        class(integrator) :: self
        procedure(RHS_X) :: F
        DOUBLE PRECISION, intent(inout) :: x(:)
        double precision, intent(in) :: tstop
        print *, 'please override advance() !'
        stop
    end subroutine advance


    subroutine cleanup(self)
        implicit none
        class(integrator) :: self

        deallocate ( self % atol)
    end subroutine cleanup

end module integrator_mod
    ! vim:et:softtabstop=4:sw=4
