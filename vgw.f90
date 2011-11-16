module vgw_mod
    use utils
    use lsode_mod
    implicit none

    type, public, abstract :: vgw
        integer :: Natom, NGAUSS
        integer :: NEQ, NG, NQNK

        character(8) :: species

        logical :: mm = .FALSE.

        double precision :: T, BL, RC, vgw_atol(3), LJA(10), LJC(10), rt
        double precision, allocatable :: mass(:), invmass(:)
        double precision, allocatable :: Y(:)

        type(lsode) :: dlsode

        integer ::  nnbmax
        integer, allocatable :: NBIDX(:,:), NNB(:)
    contains
        procedure(vgw_initX), deferred :: init
        procedure :: Ueff => vgw_Ueff
        procedure(vgw_init_Qnk0X), deferred :: init_Qnk0
        procedure(vgw_propagateX), deferred :: propagate
        procedure(vgw_logdetX), deferred :: logdet
        procedure(vgw_cleanupX), deferred :: cleanup
        procedure :: set_bl => vgw_set_bl
        procedure :: get_q => vgw_get_q
        procedure :: set_rc => vgw_set_rc
        procedure :: set_mass => vgw_set_mass
        procedure :: interaction_lists => vgw_interaction_lists
        procedure :: set_mm => vgw_set_mm
        procedure :: common_init => vgw_common_init
        procedure :: common_destroy => vgw_common_destroy
        procedure , nopass :: get_Qnk_ptr => vgw_get_Qnk_ptr
    end type vgw

    abstract interface
        subroutine vgw_initX(this, Natom, species)
            import vgw
            class(vgw) :: this
            integer, intent(in) :: Natom
            character(*), intent(in) :: species
        end subroutine

        subroutine vgw_cleanupX(this)
            import vgw
            class(vgw) :: this
        end subroutine

        subroutine vgw_init_Qnk0X(this)
            import vgw
            class(vgw) :: this
        end subroutine

        double precision function vgw_logdetX(this)
            import vgw
            class(vgw) :: this
        end function

        subroutine vgw_propagateX(this, tau)
            import vgw
            class(vgw) :: this
            double precision, intent(in) :: tau
        end subroutine
    end interface
contains
    subroutine vgw_common_init(this)
        class(vgw) :: this

        this%NEQ = 3 * this%Natom + this%NG + 1

        allocate ( this%NNB( this%Natom ), this%NBIDX(this%Natom,this%Natom), &
            this%y ( this%NEQ ), &
            this%mass(this%Natom), this%invmass( this%Natom ) )

        this%vgw_atol = 1d-4
        call this%dlsode%set_dt(0d0, 0d0, 0d0)

        if (this%species=='pH2-4g') then
            this%NGAUSS=4
            this%LJA(1:4) = (/ 1.0382522151, 0.59740391094, 0.1964765722778, &
                        0.066686117717 /)
            this%LJC(1:4) = (/ 96609.4882898, 14584.620755075, -365.4606149565, &
                        -19.55346978000 /)
            this%mass = 2.0*0.020614788876D0
            this%rc = 8.0

            this%vgw_atol = (/ 1d-3, 1d-4, 1d0 /)

            call this%dlsode%set_dt(5d-4, 1d-5, 2d-3)
        else if (this%species == 'LJ') then
            this%NGAUSS = 3
            this%LJA(1:3) = (/ 6.65, 0.79, 2.6 /)
            this%LJC(1:3) = (/ 1840d0, -1.48d0, -23.2d0 /)
            this%mass = 1.0
            this%rc = 2.5
            !this%rfullmatsq = 1.8d0**2

            this%vgw_atol = (/ 1d-5, 1d-7, .1d0 /)

            call this%dlsode%set_dt(1d-5, 1d-7, 0.25d0)
        else if (this%species == 'OH') then
            this%NGAUSS=5
            this%LJA(1:5) = (/-0.0512, 0.7628, 3.5831, 16.4327, 0.0/)
            this%LJC(1:5) = (/2.1352, 0.2919, 0.1067, 0.0571, -2.3864/)
            this%rc = 100d0

            this%mass(1) = 16d0*1823d0
            this%mass(2) = 1d0*1823d0

            if (size(this%mass) > 2) then
                write(*,*) 'Error: can only do OH dimers!'
            end if

            this%vgw_atol = (/ 1d-5, 1d-7, .1d0 /)

            call this%dlsode%set_dt(1d-5, 1d-7, 0.25d0)
            call this%dlsode%set_dt(1d-5, 1d-7, 10d0)
        end if
        this%invmass = 1d0/this%mass
    end subroutine


    subroutine vgw_common_destroy(this)
        class(vgw) :: this

        deallocate(this%y, this%mass, this%invmass, this%nnb, this%nbidx)
        call this%dlsode%destroy()
    end subroutine


    subroutine vgw_interaction_lists(this, Q)
        implicit none
        class(vgw) :: this
        real*8, intent(in) :: Q(:,:)
        integer :: N,I,J, NN
        real*8 rsq,rc2,qij(3)

        N = size(Q, 2)
        rc2=this%rc**2

        this%NBIDX = 0
        this%NNB = 0
        do I=1,N-1
            NN = 0
            do J=I+1,N
                qij=Q(:,I)-Q(:,J)
                rsq = sum(min_image(qij, this%BL)**2)
                if(rsq <= rc2) then
                    NN = NN + 1
                    this%NBIDX(NN, I) = J
                endif
            enddo
            this%NNB(i) = NN
        enddo

        this%nnbmax = maxval(this%nnb)
    end subroutine

    SUBROUTINE vgw_Ueff(this, Q0, beta, Ueff, WX)
        IMPLICIT NONE
        class(vgw) :: this
        double precision, intent(in) :: Q0(:,:), beta
        double precision, intent(out), optional :: Ueff, WX(:,:)
        real*8 :: logrho, TSTOP, start_time, stop_time
        integer :: ncalls

        call this%interaction_lists(Q0)


        this%y  = this%vgw_atol(2)
        ! set up the absolute tolerance
        this%y (1 : 3 * this%Natom) = this%vgw_atol(1)
        this%y (3 * this%Natom + 1 : 3 * this%Natom + this%NG) = this%vgw_atol(2)
        ! tolerance for Qnkp and gamakp
        this%y( this%NEQ ) = this%vgw_atol(3)

        call this%dlsode%init( this%NEQ )
        call this%dlsode%set_atol(this%y)
        
        ! initialize q_0 and Qnk
        this%y = 0d0
        this%y(1 : 3 * this%Natom) = reshape(Q0, (/ 3 * this%Natom /) )
        if (this%mm) then
            call this%init_Qnk0()
        end if

        ! solve the VGW equations, measure CPU time
        call cpu_time(start_time)

        this%T = 0
        TSTOP = 0.5d0*beta
        
        call this%propagate(0.5d0 * beta)

        call cpu_time(stop_time)

        ! write statistics
        !write (*,*) this%dlsode%get_nsteps(), 'steps,', &
        !    this%dlsode%get_ncalls(), ' RHSS calls'

        if (present(Ueff) .and. this%T > 0) then
            logrho = 2.0 * this%Natom * this%y( this%NEQ ) &
                - 0.5*this%logdet() - 1.5 * this%Natom * log(4.0*M_PI)
        
            Ueff = -logrho/beta
        end if

        if (present(WX)) then
            WX = -2d0/beta &
                * reshape(this%y(this%NEQ - 3 * this%Natom : this%NEQ - 1), &
                (/3, this%Natom/) )
        end if

        this%rt = 0
        if ( this%dlsode%get_ncalls() > 0) then
            this%rt = (stop_time - start_time) / real(this%dlsode%get_ncalls())
        end if
    end subroutine

    function vgw_classical_Utot(this, Q0) result(U)
        implicit none
        class(vgw) :: this
        double precision, intent(in) :: Q0(:,:)
        double precision :: U
        INTEGER  I,J,N
        real*8 :: rsq, QIJ(3)

        N = size(Q0, 2)

        U=0d0
        DO I=1,N-1
            DO J=I+1,N
                    qij = Q0(:,I) - Q0(:,J)
                    rsq = sum(min_image(qij, this%BL)**2)
                    U = U + sum( this%LJC (1 : this%NGAUSS) * &
                        EXP(-this%LJA ( 1 : this%NGAUSS ) * rsq) )
            ENDDO
        ENDDO
    end function


    subroutine vgw_set_bl(this, bl)
        class(vgw) :: this
        double precision :: bl

        this%bl = bl
    end subroutine

    function vgw_get_q(this)
        implicit none
        class(vgw) :: this
        double precision :: vgw_get_q(3, this%Natom)

        vgw_get_q = reshape(this%y ( 1 : 3 * this%Natom), (/ 3, this%Natom /) )
    end function


    subroutine vgw_set_rc(this, rc)
        implicit none
        class(vgw) :: this
        double precision :: rc

        this%rc = rc
    end subroutine


    subroutine vgw_set_mass(this, m)
        implicit none
        class(vgw) :: this
        double precision, intent(in) :: m(:)

        this%mass = m
        this%invmass = 1d0/m
    end subroutine


    subroutine vgw_set_mm (this, mm)
        implicit none
        class(vgw) :: this
        logical :: mm

        this%mm = mm

        this%NEQ = 3 * this%Natom + this%NG + 1
        if (this%mm) then
            this%NEQ = this%NEQ + this%NQnk + 3*this%Natom
        end if

        if (size(this%y) /= this%NEQ) then
            deallocate(this%y)
            allocate( this%y( this%NEQ ))
        end if

    end subroutine

    subroutine vgw_get_Qnk_ptr(v, ptr)
        class(vgw), pointer :: v
        double precision, pointer :: ptr(:)

        integer :: skip
        
        skip = 3*v%Natom + v%NG
        ptr => v%y(skip + 1 : skip + v%NQnk)
    end subroutine


end module
