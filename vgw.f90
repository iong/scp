module vgw_mod
    use utils
    use lsode_mod
    implicit none

    type, public, abstract :: vgw
        integer :: Natom, NGAUSS
        integer :: NEQMAX, NG, NQNK

        double precision :: T, BL, RC, vgw_atol(3), LJA(10), LJC(10), rt
        double precision, allocatable :: Y(:), mass(:), invmass(:)

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
        procedure :: set_species => vgw_set_species
        procedure :: set_rc => vgw_set_rc
        procedure :: set_mass => vgw_set_mass
        procedure :: interaction_lists => vgw_interaction_lists
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

        subroutine vgw_propagateX(this, T)
            import vgw
            class(vgw) :: this
            double precision, intent(in) :: T
        end subroutine
    end interface
contains
    subroutine vgw_set_bl(this, bl)
        class(vgw) :: this
        double precision :: bl

        this%bl = bl
    end subroutine

    subroutine vgw_get_q(this, q)
        implicit none
        class(vgw) :: this
        double precision, intent(out) :: q(:,:)

        q = reshape(this%y ( 1 : 3 * this%Natom), (/ 3, this%Natom /) )
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

    SUBROUTINE vgw_Ueff(this, Q0, beta, Ueff, WX)
        IMPLICIT NONE
        class(vgw) :: this
        double precision, intent(in) :: Q0(:,:), beta
        double precision, intent(out) :: Ueff
        double precision, intent(out), optional :: WX(:,:)
        real*8 :: logrho, TSTOP, start_time, stop_time
        integer :: i, j, ncalls, info, NEQ

        call interaction_lists(Q0)

        NEQ = 3 * this%Natom + this%NG + 1

        if (present(WX)) then
            NEQ = NEQ + this%NQnk + 3*this%Natom
        end if

        if (size(this%y) /= NEQ) then
            deallocate(this%y)
            allocate( this%y( NEQ ))
        end if


        ! set up the absolute tolerance
        this%y (1 : 3 * this%Natom) = this%vgw_atol(1)
        this%y (3 * this%Natom + 1 : 3 * this%Natom + this%NG) = this%vgw_atol(2)
        ! tolerance for Qnkp and gamakp
        this%y(NEQ) = this%vgw_atol(3)

        call this%dlsode%init(NEQ)
        call this%dlsode%set_atol(this%y)
        
        ! initialize q_0 and Qnk
        this%y = 0d0
        this%y(1 : 3 * this%Natom) = reshape(Q0, (/ 3 * this%Natom /) )
        if (present(WX)) then
            call this%init_Qnk0()
        end if

        ! solve the VGW equations, measure CPU time
        call cpu_time(start_time)

        this%T = 0
        TSTOP = 0.5d0*beta
        
        call this%propagate(0.5d0 * beta)

        logrho = 2.0 * this%Natom * this%y(NEQ) - 0.5*this%logdet() &
            - 1.5 * this%Natom * log(4.0*M_PI)
        
        call cpu_time(stop_time)

        ! write statistics
        write (*,*) this%dlsode%get_nsteps(), 'steps,', &
            this%dlsode%get_ncalls(), ' RHSS calls'

        Ueff = -logrho/beta

        if (present(WX)) then
            WX = -2d0/beta *reshape(this%y(NEQ - 3 * this%Natom : NEQ - 1), &
                (/3, this%Natom/) )
        end if

        this%rt = (stop_time - start_time) / real(ncalls)

        call this%dlsode%destroy()
    end subroutine


    subroutine vgw_set_species(this, species)
        class(vgw) :: this
        character(*) :: species

        this%vgw_atol = 1d-4
        call this%dlsode%set_dt(0d0, 0d0, 0d0)

        if (species=='pH2-4g') then
            this%NGAUSS=4
            this%LJA(1:4) = (/ 1.0382522151, 0.59740391094, 0.1964765722778, &
                        0.066686117717 /)
            this%LJC(1:4) = (/ 96609.4882898, 14584.620755075, -365.4606149565, &
                        -19.55346978000 /)
            this%mass = 2.0*0.020614788876D0
            this%rc = 8.0

            this%vgw_atol = (/ 1d-3, 1d-4, 1d0 /)

            call this%dlsode%set_dt(5d-4, 1d-5, 2d-3)
        else if (species == 'LJ') then
            this%NGAUSS = 3
            this%LJA(1:3) = (/ 6.65, 0.79, 2.6 /)
            this%LJC(1:3) = (/ 1840d0, -1.48d0, -23.2d0 /)
            this%mass = 1.0
            this%rc = 2.5
            !this%rfullmatsq = 1.8d0**2

            this%vgw_atol = (/ 1d-5, 1d-7, .1d0 /)

            call this%dlsode%set_dt(1d-5, 1d-7, 0.25d0)
        else if (species == 'OH') then
            this%NGAUSS=5
            this%LJA(1:5) = (/-0.0512, 0.7628, 3.5831, 16.4327, 0.0/)
            this%LJC(1:5) = (/2.1352, 0.2919, 0.1067, 0.0571, -2.3864/)
            this%rc = 100d0

            this%mass(1) = 16d0
            this%mass(2) = 1d0

            if (size(this%mass) > 2) then
                write(*,*) 'Error: can only do OH dimers!'
            end if

            this%vgw_atol = (/ 1d-5, 1d-7, .1d0 /)

            call this%dlsode%set_dt(1d-5, 1d-7, 0.25d0)
        end if
    end subroutine

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


    subroutine vgw_interaction_lists(this, Q)
        implicit none
        class(vgw) :: this
        real*8, intent(in) :: Q(:,:)
        integer :: N,I,J, NN
        real*8 rsq,rc2,qij(3)

        N = size(Q, 2)
        rc2=this%rc**2

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

end module
