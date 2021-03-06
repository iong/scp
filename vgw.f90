module vgw_mod
    use utils
!    use integration
    implicit none
    private

    integer, parameter :: info = 1
    integer :: debug = 0

    type, public :: vgw


        integer :: Natom, N3, NGAUSS
        integer :: NEQ, NG, NQNK, niterations

        character(8) :: species

        logical :: mm = .FALSE.

        double precision :: U, kT, BL, LJA(10), LJC(10), rt,  sigma0, epsilon0
        double precision :: q_atol, g_atol, dt_min, dt0, dt_max, qconv, gconv
        double precision :: iter_time, logdet_time
        double precision, allocatable :: Y(:), UX(:)

        !class(integrator), pointer :: prop
    contains
        procedure :: init
        procedure :: cleanup
        procedure :: F => vgw_F
        procedure :: analyze
        procedure :: init_prop
        procedure :: logdet
        procedure :: converge
        procedure :: regtransrot
        procedure :: rhs
        procedure :: set_bl
        procedure :: get_q
        procedure :: Utot0
    end type vgw

    public :: transrot_subspace, info, debug

 contains
    subroutine init(self, Natom, species)
        class(vgw) :: self
        integer, intent(in) :: Natom
        character(*), intent(in) :: species
        double precision, parameter :: xeqOH = 1.8897, kOH = 0.49536d0

        self%species = species
        self%Natom = Natom
        self%N3 = 3 * Natom
        self%niterations = 0
        self % iter_time = 0d0
        self % logdet_time = 0d0

        allocate (self % UX (3*self % Natom))

        self % q_atol = 1d-4
        self % g_atol = 1d-4
        self % dt_min = 0d0
        self % dt_max = 0d0
        self % dt0    = 0d0

        self % sigma0 = 1d0
        self % epsilon0 = 1d0
        self % bl = 1d10

        if (self % species=='pH2-4g') then
            self % NGAUSS=4
            self % LJA(1:4) = (/ 1.0382522151, 0.59740391094, 0.1964765722778, &
                    0.066686117717 /)
            self % LJC(1:4) = (/ 96609.4882898, 14584.620755075, -365.4606149565, &
                    -19.55346978000 /)

            self % q_atol = 1d-3
            self % g_atol = 1d-4

            self % dt_min = 1d-5
            self % dt_max = 2d-3
            self % dt0 = 5d-4
        else if (self % species == 'LJ') then
            self % NGAUSS = 3
            self % LJA(1:3) = (/ 6.65, 0.79, 2.6 /)
            self % LJC(1:3) = (/ 1840d0, -1.48d0, -23.2d0 /)

            self % q_atol = 1d-3
            self % g_atol = 1d-5
        else if (self % species == 'OH') then
            self % NGAUSS=5
            self % LJA(1:5) = (/-0.0902, 0.9454, 4.0565, 17.5229, 0.0/)
            self % LJC(1:5) = (/2.0233, 0.5054, 0.1912, 0.1085, -2.4140/)
            self % LJA(1:5) =  self % LJA(1:5) / xeqOH**2
            self % LJC(1:5) = self % LJC(1:5) * kOH* xeqOH**2

            self % q_atol = 1d-5
            self % g_atol = 1d-7

            self % dt_min = 1d-7
            self % dt_max = 10d0
            self % dt0 = 1d-5
        else if (self % species == 'HarmOsc') then
            self % NGAUSS=5
            self % LJA(1:5) = (/-0.0854686, 0.227202, 2.17897, 8.00009,  0.0/)
            self % LJC(1:5) = (/3.42986, -1.05074, 0.0294075, -0.00688583,-2.4021/)
            self % LJA(1:5) =  self % LJA(1:5) / xeqOH**2
            self % LJC(1:5) = self % LJC(1:5) * kOH * xeqOH**2


            self % q_atol = 1d-5
            self % g_atol = 1d-7

            self % dt_min = 1d-7
            self % dt_max = 10d0
            self % dt0 = 1d-5
        end if
    end subroutine init


    subroutine cleanup(self)
        class(vgw) :: self

        !call self % prop % cleanup()
        deallocate(self % y, self % UX)!, self % prop)
    end subroutine cleanup


    subroutine analyze(self, q0)
        class(vgw), target :: self
        double precision, intent(in) :: q0(:)
    end subroutine analyze


    subroutine init_prop(self, q0)
        class(vgw), target :: self
        double precision, intent(in) :: q0(:)

        integer :: N3

        N3 = 3 * self % Natom

        if (allocated(self % y)) then
                deallocate(self%y)
        end if

        allocate(self % y ( self % NEQ ))

        self % y(1 : N3) = q0
        self % y(N3+1:) = 0d0
    end subroutine init_prop

    function logdet(self)
        class(vgw), target :: self
        double precision :: logdet

        print *, 'please override logdet!'
        logdet = 0d0
    end function logdet

    subroutine regtransrot(self, y)
        class(vgw) :: self
        double precision, intent(inout) :: y(:)
    end subroutine regtransrot

    function test_convergence(y0, y1, rtol, atol) 
        implicit none
        double precision, intent(in) :: y0(:), y1(:), rtol, atol
        logical :: test_convergence

        test_convergence = all(abs(y0-y1) < (0.5*rtol*( abs(y0) + abs(y1) ) + atol))
    end function


    subroutine converge(self, rtol, atol)
        IMPLICIT NONE
        class(vgw) :: self
        double precision, intent(in) :: rtol, atol

        double precision, allocatable :: yproj(:), yprojold(:)

        integer :: n, nmax, N3
        double precision ::  dt, t

        logical :: converged

        N3 = 3 * self%Natom
        allocate(yproj(self%NEQ), yprojold(self%NEQ))

        dt = 0.0125d0
        t = 0d0
        yproj = 0d0
        do
            yprojold = yproj
            self%y = self%y + dt * self % rhs(self%y, t)
            
           ! call self%regtransrot(self % y)
            yproj = self % y
            converged = test_convergence(yproj(:N3), yprojold(:N3), 0d0, dt*1d-2) &
                    .AND. test_convergence(yproj(N3+1:), yprojold(N3+1:), rtol, dt*atol)

            t = t + dt
            self % niterations = self % niterations + 1

            ! reduce the time step after the two orders of magnitude of the
            ! error have been eliminated.
            !if (err < 0.01d0/rtol) dt = dt / 2d0

            if (converged) exit
        end do

        deallocate(yproj, yprojold)
    end subroutine CONVERGE


    function rhs(self, Y, T)
        IMPLICIT NONE
        class(vgw) :: self
        double precision, intent(in) :: T
        double precision, intent(in), target :: Y(:)
        double precision, target :: rhs(size(Y))
    end function



    
    function vgw_F(self, Q0, kT, noinit)
        IMPLICIT NONE
        class(vgw) :: self
        double precision, intent(in) :: Q0(:), kT
        logical, intent(in), optional :: noinit
        double precision :: vgw_F
        real*8 ::  t1, t2, t3, c0
        integer :: i
        logical :: lnoinit = .FALSE.

        self % kT = kT / self % epsilon0
        c0 = 2d0 * M_PI * exp(0.5) *self%epsilon0 * self%sigma0

        if (present(noinit)) lnoinit = noinit
        if (.NOT. lnoinit) then
            call self % analyze(q0 / self % sigma0)
            call self % init_prop(q0 / self % sigma0)
        end if

        ! solve the VGW equations, measure CPU time
        call cpu_time(t1)
        call self % converge(1d-4, 1d-5)
        call cpu_time(t2)
        vgw_F =  (self%U - self%kT * self%logdet() )/self%Natom &
                    - 3  * self%kT * log(self % kT * c0)
        vgw_F = self%epsilon0 * vgw_F
        call cpu_time(t3)

        self % iter_time = self % iter_time + t2 - t1
        self % logdet_time = self % logdet_time + t3 - t2
    end function vgw_F


    function Utot0(self, Q0) result(U)
        implicit none
        class(vgw) :: self
        double precision, intent(in) :: Q0(:,:)
        double precision :: U
        INTEGER  I,J,N
        real*8 :: rsq, QIJ(3)

        N = size(Q0, 2)

        U=0d0
        DO I=1,N-1
            DO J=I+1,N
                qij = Q0(:,I) - Q0(:,J)
                rsq = sum(min_image(qij, self % BL)**2)
                U = U + sum( self % LJC (1 : self % NGAUSS) * &
                        EXP(-self % LJA ( 1 : self % NGAUSS ) * rsq) )
            ENDDO
        ENDDO
    end function Utot0


    subroutine set_bl(self, bl)
        class(vgw) :: self
        double precision :: bl

        self % bl = bl
    end subroutine set_bl


    function get_q(self)
        implicit none
        class(vgw) :: self
        double precision :: get_q(3, self % Natom)

        get_q = reshape(self % y ( 1 : 3 * self % Natom), (/ 3, self % Natom /) )
    end function get_q

    function transrot_subspace(q) result(U)
        implicit none
        double precision, intent(in) :: q(:)
        double precision :: U(size(q),6)

        integer :: i, j, N

        N = size(q)

        U = 0d0
        U(1::3,1) = 1d0/sqrt(N/3d0)
        U(2::3,2) = 1d0/sqrt(N/3d0)
        U(3::3,3) = 1d0/sqrt(N/3d0)
        U(2::3,4) = -q(3::3)
        U(3::3,4) =  q(2::3)
        U(1::3,5) =  q(3::3)
        U(3::3,5) = -q(1::3)
        U(1::3,6) = -q(2::3)
        U(2::3,6) =  q(1::3)

        do i=4,6
            do j=1,i-1
                U(:,i) = U(:,i) - U(:,j) * sum ( U(:,i) * U(:,j) )
            end do
            U(:,i) = U(:,i) / sqrt(sum(U(:,i)**2))
        end do
    end function transrot_subspace
end module vgw_mod
