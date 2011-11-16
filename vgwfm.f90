module vgwfm_mod
    use vgw_mod
    implicit none

    type, public, extends(vgw) :: vgwfm
    real*8 :: rfullmatsq

    real*8, allocatable :: G(:,:), GP(:,:), GU(:,:)
    real*8, allocatable :: UX(:), UXY(:,:)
contains
    procedure :: init => vgwfm_init
    procedure :: cleanup => vgwfm_cleanup
    procedure :: init_Qnk0 => vgwfm_init_Qnk0
    procedure :: propagate => vgwfm_propagate
    procedure :: logdet => vgwfm_logdet
    procedure :: get_Meff => vgwfm_get_Meff
end type vgwfm

contains

subroutine vgwfm_init(this, Natom, species)
    implicit none
    class(vgwfm) :: this
    integer, intent(in) :: Natom
    character(*), intent(in) :: species

    this%Natom = Natom
    this%NG =  (9*Natom**2 + 3*Natom)/2
    this%NQNK = 9*Natom**2

    this%species = species

    allocate(this%G(3*Natom,3*Natom), this%GP(3*Natom, 3*Natom), &
            this%GU(3*Natom, 3*Natom), this%UX(3*Natom), &
            this%UXY(3*Natom, 3*Natom))

    call this%common_init()
end subroutine vgwfm_init

subroutine vgwfm_cleanup(this)
    implicit none
    class(vgwfm) :: this

    deallocate(this%G, this%GP, this%GU, this%UX, this%UXY)

    call this%common_destroy()
end subroutine vgwfm_cleanup

subroutine fm_get_g(N, y, G)
    implicit none
    integer, intent(in) :: N
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: G(:,:)

    integer :: ypos, i

    ypos = 3 * N + 1
    do i = 1, 3 * N
        G( i : 3 * N, i) = y (ypos : ypos + 3*N - i)
        G( i, i : 3 * N) = y (ypos : ypos + 3*N - i)
        ypos = ypos + 3 * N - i + 1
    end do
end subroutine fm_get_g

subroutine fm_set_g(N, G, y)
    implicit none
    integer, intent(in) :: N
    double precision, intent(in) :: G(:,:)
    double precision, intent(out) :: y(:)

    integer :: ypos, i

    ypos = 3 * N + 1
    do i = 1,3 * N
        y (ypos : ypos + 3 * N - i) = G(i : 3 * N, i)
        ypos = ypos + 3 * N - i + 1
    end do
end subroutine fm_set_g


subroutine vgwfm_init_Qnk0(this)
    class(vgwfm) :: this
    integer :: j, skip, jdiag

    skip = 3 * this%Natom + this%NG
    jdiag = 1
    do j=1, 3*this%Natom
        this%y (skip + jdiag) = 1d0
        jdiag = jdiag + 3*this%Natom + 1
    end do
end subroutine vgwfm_init_Qnk0


function vgwfm_logdet(this)
    class(vgwfm) :: this
    double precision :: vgwfm_logdet

    integer :: info, j

    call fm_get_g(this%Natom, this%y, this%G)
    call dpotrf('U', 3 * this%Natom, this%G, 3 * this%Natom, info)

    vgwfm_logdet=0.0
    DO j=1, 3 * this%Natom
        vgwfm_logdet = vgwfm_logdet + LOG(ABS( this%G(j,j) ))
    ENDDO
    vgwfm_logdet = 2d0* vgwfm_logdet
end function vgwfm_logdet


subroutine vgwfm_propagate(this, tau)
    class(vgwfm) :: this
    double precision :: tau

    call this%dlsode%propagate(RHSSFM, this%y, this%T, tau)

contains

    SUBROUTINE RHSSFM(NEQ, T, Y, YP)
        !    use omp_lib
        IMPLICIT NONE
        integer, intent(in) :: NEQ!, IPAR(:)
        double precision, intent(in) :: T
        double precision, intent(in), target :: Y(NEQ)!, RPAR(:)
        double precision, intent(out), target :: YP(NEQ)
        INTEGER :: J,I1,I2,IG, I1_3, I2_3, skip
        REAL*8 AG(3,3), &
                DETA,DETAG,QZQ,U12, &
                G12(3,3),A(3,3), &
                Zq(3), Z(3,3),Q12(3)
        real*8 :: UXY0(3,3), UX0(3), U, TRUXXG, k1

        double precision, pointer :: Qnk(:), Qnkp(:), gamak(:)

        if (y(3 * this%Natom + 1) == 0d0) then
            call rhss_zero_time(NEQ, y, yp)
            return
        end if

        !    do I1=1,Natom-1
        call fm_get_g(this%Natom, this%y, this%G) 

        U = 0; this%UX = 0; this%UXY = 0;
        do I1 = 1, this%Natom - 1
            I1_3 = 3*(I1-1) + 1
            do I2 = I1+1, this%Natom
                I2_3 = 3*(I2-1) + 1
                Q12 = y(I1_3:I1_3+2) - y(I2_3:I2_3+2)
                Q12 = min_image(Q12, this%bl)
                G12=this%G (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) &
                        + this%G (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) &
                        - this%G (I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) &
                        - this%G (I1_3 : I1_3 + 2, I2_3 : I2_3 + 2)

                call detminvm(G12, DETA, A)
                DETA = 1.0d0/DETA

                UX0 = 0d0; UXY0 = 0d0
                DO IG=1, this%NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
                    AG = A
                    do J=1,3
                        AG(J,J) = this%LJA(IG) + AG(J,J)
                    end do

                    call detminvm(AG, DETAG, Z)
                    Z = - this%LJA(IG)**2 * Z

                    do J=1,3
                        Z(J,J) = Z(J,J) + this%LJA(IG)
                    end do

                    Zq = matmul(Z, Q12) ! R = -2.0*Zq
                    qZq = dot_product(Q12, Zq) 

                    U12 = SQRT(DETA/DETAG) * EXP(-qZq) * this%LJC(IG)
                    U = U + U12

                    UX0 = UX0 - 2d0*U12*Zq
                    do J=1,3
                        UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                    end do
                end do ! IG

                ! Avoid load and store as much as possbile. Store now, process as a stream later. Much faster.
                this%UX (I1_3 : I1_3 + 2) = this%UX (I1_3 : I1_3 + 2) + UX0
                this%UX (I2_3 : I2_3 + 2) = this%UX (I2_3 : I2_3 + 2) - UX0

                this%UXY (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) = &
                        this%UXY (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) + UXY0

                this%UXY (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) = &
                        this%UXY (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) + UXY0

                !if (sum(q12**2) <= rfullmatsq) then
                this%UXY (I1_3 : I1_3 + 2, I2_3 : I2_3 + 2) = -UXY0
                this%UXY (I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) = -UXY0
                !end if
            end do ! I2
        end do ! I1


        !QP = -matmul(G, UX)
        call dsymv('L', 3 * this%Natom, -1d0, this%G, 3 * this%Natom, this%UX, 1, &
                0d0, yp, 1)

        TRUXXG = sum( this%UXY * this%G)

        !GU = matmul(G, UXY)
        call dsymm('L', 'L', 3 * this%Natom, 3 * this%Natom, 1d0, this%G, 3 * this%Natom, &
                this%UXY, 3 * this%Natom, 0d0, this%GU, 3 * this%Natom)

        !GP = -matmul(GU, G)
        call dsymm('R', 'L', 3 * this%Natom, 3 * this%Natom, -1d0, this%G, 3 * this%Natom,&
                this%GU, 3 * this%Natom, 0d0, this%GP, 3 * this%Natom)

        do i1=1, this%Natom
            k1=             3*(i1 - 1) + 1

            this%GP (k1, k1)     = this%GP (k1, k1)     + this%invmass(i1)
            this%GP (k1+1, k1+1) = this%GP (k1+1, k1+1) + this%invmass(i1)
            this%GP (k1+2, k1+2) = this%GP (k1+2, k1+2) + this%invmass(i1)
        end do

        call fm_set_g(this%Natom, this%GP, yp)

        if (this%mm) then
            skip = 3*this%Natom + this%NG
            Qnk => y( skip + 1 : skip + this%NQnk )
            Qnkp => yp( skip + 1 : skip + this%NQnk )

            call dgemm('N', 'N', 3*this%Natom, 3*this%Natom, 3*this%Natom, -1d0, this%GU, &
                    3*this%Natom, Qnk, 3*this%Natom, 0d0, Qnkp, 3*this%Natom)

            skip = skip + this%NQnk
            gamak => yp(skip + 1 : skip + 3*this%Natom)
            call dgemv('T', 3 * this%Natom, 3 * this%Natom, -1d0, Qnk, 3 * this%Natom, &
                    this%UX, 1, 0d0, gamak, 1)
        end if

        yp(NEQ) = -(0.25d0*TRUXXG + U)/real(this%Natom)
    END SUBROUTINE RHSSFM


    subroutine rhss_zero_time(NEQ, y, yp)
        integer, intent(in) :: NEQ
        double precision, intent(in) :: y(:)
        double precision, intent(out) :: yp(:)

        double precision :: qij(3), qi(3), qj(3), rsq, LJC_EXPAV(this%NGAUSS), &
                UX0(3), U
        integer :: i, j, ki, kj

        yp = 0d0

        j = 3 * this%Natom + 1
        do i=1, 3 * this%Natom
            yp(j) = this%invmass( (i-1)/3 + 1)
            j = j + 3 * this%Natom - i + 1
        end do

        U=0d0

        DO I=1, this%Natom-1
            qi = y(3*I-2:3*I)
            DO J=1, this%NNB(I)
                qj = y(3 * this%NBIDX(J,I)-2 : 3 * this%NBIDX(J,I))
                qij = qi - qj
                qij = min_image(qij, this%BL)
                rsq = sum(qij**2)

                LJC_EXPAV = this%LJC (1 : this%NGAUSS) * EXP(-rsq &
                        * this%LJA (1 : this%NGAUSS))

                U = U + sum(LJC_EXPAV)

                ! there are to -1 here which cancel each other
                UX0 = qij*(2.0*sum( this%LJA (1 : this%NGAUSS) * LJC_EXPAV ) )


                ki = 3*this%Natom + this%NG + this%NQnk + 3*I  - 2
                kj = 3*this%Natom + this%NG + this%NQnk + 3*this%NBIDX(J,I)  - 2

                yp(ki : ki + 2) = yp(ki : ki + 2) + UX0
                yp(kj : kj + 2) = yp(kj : kj + 2) - UX0
            ENDDO
        ENDDO

        yp(NEQ) = -U/real(this%Natom)
    end subroutine rhss_zero_time
end subroutine vgwfm_propagate

subroutine vgwfm_get_Meff(this, Meff, invMeff, sqrtMeff, sqrtInvMeff)
    class(vgwfm) :: this
    double precision, intent(out) :: Meff(:,:)
    double precision, intent(out), optional :: invMeff(:,:), sqrtMeff(:,:), &
            sqrtInvMeff(:,:)

    double precision, dimension(3 * this%Natom) :: W, FV1, FV2
    double precision, dimension(3 * this%Natom, 3*this%Natom) :: U, UD
    integer :: j, ierr, N

    N = 3 * this%Natom

    call fm_get_g(this%Natom, this%y, Meff)

    do j=1,N
        W(j) = this%mass(1 + (j-1)/3)
    end do

    ! M*G*M
    do j=1,N
        Meff(:,j) = this%T * Meff(:,j) * W * W(j)/this%T
    end do

    if (present(invMeff)) then
        call RS(N,N,Meff,W,1,U,FV1,FV2,ierr)
        do j=1,N
            UD(:,j) = (1d0/W(j))*U(:,j)
        end do
        call dgemm('N', 'T', N, N, N, 1d0, U, N, UD, N,  0d0, invMeff, N)

        if (present(sqrtMeff)) then
            do j=1,N
                UD(:,j) = sqrt(W(j))*U(:,j)
            end do
            call dgemm('N', 'T', N, N, N, 1d0, U, N, UD, N,  0d0, sqrtMeff, N)
        end if

        if (present(sqrtInvMeff)) then
            do j=1,N
                UD(:,j) = sqrt(1d0/W(j))*U(:,j)
            end do
            call dgemm('N', 'T', N, N, N, 1d0, U, N, UD, N,  0d0, sqrtInvMeff, N)
        end if
    end if
end subroutine vgwfm_get_Meff
end module vgwfm_mod
