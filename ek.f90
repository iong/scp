module ek_mod
    use integrator_mod
    private

    type, public, extends(integrator) ::  ek
        double precision, allocatable :: xe(:), x0(:), x1(:), xp(:,:)
    contains
        procedure :: init
        procedure :: advance
        procedure :: step 
        procedure :: cleanup
    end type ek

contains
    subroutine  init(self,F,  NEQ, tstart, dt0)
        implicit none
        class(ek) :: self
        procedure(RHS_X) :: F
        integer, intent(IN) :: NEQ
        double precision, intent(in) :: tstart, dt0

        call self % integrator % init(F, NEQ, tstart, dt0)

        allocate(self%xe(NEQ), self%x0(NEQ), self%x1(NEQ), self%xp(NEQ, 2))
    end subroutine init

    subroutine cleanup(self)
        implicit none
        class(ek) :: self
        deallocate(self%xe, self%x0, self%x1, self%xp)
        call self % integrator % cleanup()
    end subroutine cleanup

    subroutine advance(self, x, tstop)
        implicit none
        class(ek) :: self
        DOUBLE PRECISION, intent(inout) :: x(:)
        double precision, intent(in) :: tstop
        DOUBLE PRECISION :: rmserr, newdt


        if (size(x) /= size(self%x1) ) then
            write (*,*) 'EK error, the size of the vector to be propagated differs', &
                    ' from internal storage.'
        end if

        self % x0 = x
        do while (self % t < tstop)
            rmserr = self % step(x, self % dt)
            write (*,*) self % t, self % dt, rmserr
            if (rmserr <= 0.5) then
                self % x0 = x
                if (self % t + self % dt > tstop) then
                    exit
                endif
                self % t = self % t + self % dt
                if (rmserr <=0.1) then
                    newdt = self % dt * 1.52799165273419
                    if (newdt <= self %dtmax) then
                        self % dt = newdt
                    end if
                endif

                self % nsteps = self % nsteps + 1
            else
                self % dt = self % dt/1.90779652734191
                x = self % x0
            endif
        enddo
        rmserr = self % step(x, tstop - self % t)
        write (*,*) self % t,  tstop - self % t, rmserr
        self % t = tstop
    end subroutine advance


    function step(self, x, dt)
        implicit none
        class(ek) :: self
        DOUBLE PRECISION, intent(inout) :: x(:)
        DOUBLE PRECISION, intent(in) :: dt
        double precision :: step

        integer, parameter :: p=2
        integer :: NEQ



        NEQ = size(x)
        call self % F(NEQ, self % t, x, self%xp(:,1))
        self % x1 = x + dt*self % xp(:,1)
        call self % F(NEQ, self % t, self % x1, self % xp(:,2))

        self % ncalls = self % ncalls + 2

        x = x + 0.5*dt*( self % xp(:,1) + self % xp(:,2) )
        self % xe = 0.5*dt*( self % xp(:,1) - self % xp(:,2) )

        self % xe = self % xe / (abs(self % rtol * x) + self % atol)
        step = sqrt( sum(self % xe **2) / NEQ )

    end function step
end module
