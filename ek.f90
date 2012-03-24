module ek_mod
    use integrator_mod
    private

    type, public, extends(integrator) ::  ek
        double precision, allocatable :: xe(:), x0(:), x1(:), xp(:,:)
    contains
        procedure :: init
        procedure :: advance
        procedure :: get_xp
        procedure :: step 
        procedure :: cleanup
    end type ek

contains
    subroutine  init(self, NEQ)
        implicit none
        class(ek) :: self
        integer, intent(IN) :: NEQ

        call self % integrator % init(NEQ)

        allocate(self%xe(NEQ), self%x0(NEQ), self%x1(NEQ), self%xp(NEQ, 2))
    end subroutine init

    subroutine cleanup(self)
        implicit none
        class(ek) :: self
        deallocate(self%xe, self%x0, self%x1, self%xp)
        call self % integrator % cleanup()
    end subroutine cleanup

    subroutine advance(self, F, x, tstop)
        implicit none
        class(ek) :: self
        procedure(RHS_X) :: F
        DOUBLE PRECISION, intent(inout) :: x(:)
        double precision, intent(in) :: tstop
        DOUBLE PRECISION :: rmserr, t0


        if (size(x) /= size(self%x1) ) then
            write (*,*) 'EK error, the size of the vector to be propagated differs', &
                    ' from internal storage.'
        end if

        t0 = self%t
        do while (self % t < tstop)
            t0 = self%t
            self % x0 = x
            rmserr = self % step(F, x, self % dt)
            if (rmserr > 0.5 .AND. self%dt > self%dtmin) then
                self % dt = max(self % dt/1.90779652734191, self % dtmin)
                x = self % x0
            elseif (rmserr <= 0.5) then
                self % t = self % t + self % dt
                self % nsteps = self % nsteps + 1

                if (rmserr <=0.1) then
                    self % dt = min(self % dt * 1.52799165273419, self %dtmax)
                endif
            endif
        enddo
        if (self % task == 1) then
            rmserr = self % step(F, x, tstop - self % t)
            x = x - (x - self%x0) * (self%t - tstop)/(self%t - t0)
            self % t = tstop
        elseif (self % task == 3) then
        end if

        self % nsteps = self % nsteps + 1
    end subroutine advance

    subroutine get_xp(self, xp)
        implicit none
        class(ek) :: self
        DOUBLE PRECISION, intent(out) :: xp(:)

        xp = 0.5d0*sum(self%xp, 2)
    end subroutine

    function step(self, F, x, dt)
        implicit none
        class(ek) :: self
        procedure(RHS_X) :: F
        DOUBLE PRECISION, intent(inout) :: x(:)
        DOUBLE PRECISION, intent(in) :: dt
        double precision :: step

        integer, parameter :: p=2
        integer :: NEQ



        NEQ = size(x)
        call F(NEQ, self % t, x, self%xp(:,1))
        self % x1 = x + dt*self % xp(:,1)
        call F(NEQ, self % t, self % x1, self % xp(:,2))

        self % ncalls = self % ncalls + 2

        x = x + 0.5*dt*( self % xp(:,1) + self % xp(:,2) )
        self % xe = 0.5*dt*( self % xp(:,1) - self % xp(:,2) )

        self % xe = self % xe / (abs(self % rtol * x) + self % atol)
        step = sqrt( sum(self % xe **2) / NEQ )

    end function step
end module
