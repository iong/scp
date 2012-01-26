module vgw_mod
    use utils
    use integration
    implicit none
    private

    type, public :: vgw
        integer :: Natom, NGAUSS
        integer :: NEQ, NG, NQNK

        character(8) :: species

        logical :: mm = .FALSE.

        double precision :: U, kT, BL, RC, LJA(10), LJC(10), rt
        double precision :: q_atol, g_atol, dt_min, dt0, dt_max
        double precision, allocatable :: mass(:), invmass(:)
        double precision, allocatable :: Y(:)

        class(integrator), pointer :: prop

        !integer ::  nnbmax
        !integer, allocatable :: NBIDX(:,:), NNB(:)
    contains
        procedure :: init
        procedure :: cleanup
        procedure :: F => vgw_F
        procedure :: init_prop
        procedure :: logdet
        procedure :: propagate
        procedure :: set_bl
        procedure :: get_q
        procedure :: set_rc
        procedure :: set_mass
    end type vgw

 contains
    subroutine init(self, Natom, species)
        class(vgw) :: self
        integer, intent(in) :: Natom
        character(*), intent(in) :: species
        double precision, parameter :: xeqOH = 1.8897, kOH = 0.49536d0

        self%species = species
        self%Natom = Natom

        allocate (self % mass(self % Natom), self % invmass( self % Natom ) )

        self % q_atol = 1d-4
        self % g_atol = 1d-4
        self % dt_min = 0d0
        self % dt_max = 0d0

        if (self % species=='pH2-4g') then
            self % NGAUSS=4
            self % LJA(1:4) = (/ 1.0382522151, 0.59740391094, 0.1964765722778, &
                    0.066686117717 /)
            self % LJC(1:4) = (/ 96609.4882898, 14584.620755075, -365.4606149565, &
                    -19.55346978000 /)
            self % mass = 2.0*0.020614788876D0
            self % rc = 8.0

            self % q_atol = 1d-3
            self % g_atol = 1d-4

            self % dt_min = 1d-5
            self % dt_max = 2d-3
            self % dt0 = 5d-4
        else if (self % species == 'LJ') then
            self % NGAUSS = 3
            self % LJA(1:3) = (/ 6.65, 0.79, 2.6 /)
            self % LJC(1:3) = (/ 1840d0, -1.48d0, -23.2d0 /)
            self % mass = 1.0
            self % rc = 2.5
            !self % rfullmatsq = 1.8d0**2

            self % q_atol = 1d-5
            self % g_atol = 1d-7

            self % dt_min = 1d-7
            self % dt_max = 0.25
            self % dt0 = 1d-5

        else if (self % species == 'OH') then
            self % NGAUSS=5
            self % LJA(1:5) = (/-0.0902, 0.9454, 4.0565, 17.5229, 0.0/)
            self % LJC(1:5) = (/2.0233, 0.5054, 0.1912, 0.1085, -2.4140/)
            self % LJA(1:5) =  self % LJA(1:5) / xeqOH**2
            self % LJC(1:5) = self % LJC(1:5) * kOH* xeqOH**2
            self % rc = 100d0

            self % mass(1) = 16d0*1823d0
            self % mass(2) = 1d0*1823d0

            if (size(self % mass) > 2) then
                write(*,*) 'Error: can only do OH dimers!'
            end if

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
            self % rc = 100d0

            self % mass(1) = 16d0*1823d0
            self % mass(2) = 1d0*1823d0

            if (size(self % mass) > 2) then
                write(*,*) 'Error: can only do OH dimers!'
            end if

            self % q_atol = 1d-5
            self % g_atol = 1d-7

            self % dt_min = 1d-7
            self % dt_max = 10d0
            self % dt0 = 1d-5
        end if
        self % invmass = 1d0/self % mass
    end subroutine init


    subroutine cleanup(self)
        class(vgw) :: self

        call self % prop % cleanup()
        deallocate(self % mass, self % invmass, self % prop)
    end subroutine cleanup


    subroutine init_prop(self, q0)
        class(vgw) :: self
        double precision, intent(in) :: q0(:,:)
    end subroutine init_prop

    function logdet(self)
        class(vgw) :: self
        double precision :: logdet

        print *, 'please override logdet!'
        logdet = 0d0
    end function logdet

    subroutine propagate(self, tstop)
        IMPLICIT NONE
        class(vgw) :: self
        double precision, intent(in) :: tstop
        
    end subroutine PROPAGATE

    function vgw_F(self, Q0, kT, tstop)
        IMPLICIT NONE
        class(vgw) :: self
        double precision, intent(in) :: Q0(:,:), kT, tstop
        double precision :: vgw_F
        real*8 ::  start_time, stop_time

        self % kT = kT

        call self % init_prop(q0)

        ! solve the VGW equations, measure CPU time
        call cpu_time(start_time)
        call self % propagate(tstop)
        call cpu_time(stop_time)

        self % rt = 0
        if ( self % prop%ncalls > 0) then
            self % rt = (stop_time - start_time) / real(self % prop%ncalls)
        end if
      
        vgw_F = -kT*( self%logdet() + 3 * self%Natom * log(2d0*kT) ) / 2d0
    end function vgw_F


    function classical_Utot(self, Q0) result(U)
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
    end function classical_Utot


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


    subroutine set_rc(self, rc)
        implicit none
        class(vgw) :: self
        double precision :: rc

        self % rc = rc
    end subroutine set_rc


    subroutine set_mass(self, m)
        implicit none
        class(vgw) :: self
        double precision, intent(in) :: m(:)

        self % mass = m
        self % invmass = 1d0/m
    end subroutine set_mass


!!$    subroutine interaction_lists(self, Q)
!!$        implicit none
!!$        class(vgw) :: self
!!$        real*8, intent(in) :: Q(:,:)
!!$        integer :: N,I,J, NN
!!$        real*8 rsq,rc2,qij(3)
!!$
!!$        N = size(Q, 2)
!!$        rc2=self % rc**2
!!$
!!$        self % NBIDX = 0
!!$        self % NNB = 0
!!$        do I=1,N-1
!!$            NN = 0
!!$            do J=I+1,N
!!$                qij=Q(:,I)-Q(:,J)
!!$                rsq = sum(min_image(qij, self % BL)**2)
!!$                if(rsq <= rc2) then
!!$                    NN = NN + 1
!!$                    self % NBIDX(NN, I) = J
!!$                endif
!!$            enddo
!!$            self % NNB(i) = NN
!!$        enddo
!!$
!!$        self % nnbmax = maxval(self % nnb)
!!$    end subroutine interaction_lists


end module vgw_mod
