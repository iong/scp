module vgwfm_mod
    use vgw_mod
    implicit none

    type, public, extends(vgw) :: vgwfm
        real*8 :: rfullmatsq

        real*8, allocatable :: G(:,:), GP(:,:), GU(:,:)
        real*8, allocatable :: UX(:), UXY(:,:)
    contains
        procedure :: init => vgwfm_init
        procedure :: Ueff => vgwfm_Ueff
        procedure :: cleanup => vgwfm_cleanup
    end type vgwfm
    
contains

subroutine vgwfm_init(Natom, species)
    use omp_lib
    implicit none
    integer, intent(in) :: Natom
    character(*), intent(in) :: species

    this%Natom = Natom
    this%NG =  (9*Natom**2 + 3*Natom)/2
    this%NQNK = 9*Natom**2
    this%NEQMAX = 3*Natom + this%NG + this%NQNK + 3*Natom + 1


    allocate(this%G(3*Natom,3*Natom), this%GP(3*Natom, 3*Natom), this%&
        this%GU(3*Natom, 3*Natom), this%UX(3*Natom), &
        this%UXY(3*Natom, 3*Natom), this%NNB(Nmax), this%NBIDX(Nmax,Nmax))
    
    call this%set_species(species)
end subroutine

subroutine vgwfm_cleanup()
    deallocate(this%G, this%GP, this%&
        this%GU, this%UX, &
        this%UXY, this%NNB, this%NBIDX)
end subroutine

subroutine fm_get_g(this, G)
    implicit none
    class(vgwfm) :: this
    double precision, intent(out) :: G(:,:)

    integer :: ypos, i

    ypos = 3 * this%Natom + 1
    do i = 1, 3 * this%Natom
        G( i : 3 * this%Natom, i) = this%y (ypos : ypos + 3*this%Natom - i)
        G( i, i : 3 * this%Natom) = this%y (ypos : ypos + 3*this%Natom - i)
        ypos = ypos + 3 * this%Natom - i + 1
    end do
end subroutine

subroutine fm_set_g(G, y)
    implicit none
    double precision, intent(in) :: G(:,:)
    double precision, intent(out) :: y(:)
    
    integer :: ypos, i

    ypos = 3 * this%Natom + 1
    do i = 1,3 * this%Natom
        y (ypos : ypos + 3 * this%Natom - i) = G(i : 3 * this%Natom, i)
        ypos = ypos + 3 * this%Natom - i + 1
    end do
end subroutine


subroutine vgwfm_init_Qnk0(this)
    class(vgwfm) :: this

    do j=1, this%Natom
        this%y (3 * this%Natom + this%NG + 9*(j - 1) + 1) = 1d0
        this%y (3 * this%Natom + this%NG + 9*(j - 1) + 5) = 1d0
        this%y (3 * this%Natom + this%NG + 9* j         ) = 1d0
    end do
end subroutine


function vgwfm_logdet(this)
    class(vgwfm) :: this
    double precision :: vgwfm_logdet

    integer :: info, j

    call this%get_g(this%G)
    call dpotrf('U', 3 * this%Natom, this%G, 3 * this%Natom, info)

    vgwfm_logdet=0.0
    DO j=1, 3 * this%Natom
        vgwfm_logdet = vgwfm_logdet + LOG(ABS( this%G(j,j) ))
    ENDDO
    vgwfm_logdet = 2d0* vgwfm_logdet
end function

subroutine vgwfm_propagate(this, tau)
    class(vgwfm) :: this
    double precision :: tau

    call this%dlsode%propagate(RHSSFM, this%y, this%T, tau)
end subroutine

SUBROUTINE RHSSFM(NEQ, T, Y, YP)
!    use omp_lib
    IMPLICIT NONE
    integer, intent(in) :: NEQ!, IPAR(:)
    double precision, intent(in) :: T, Y(NEQ)!, RPAR(:)
    double precision, intent(out) :: YP(NEQ)
    INTEGER :: J,I1,I2,IG, I1_3, I2_3
    REAL*8 AG(3,3), &
            DETA,DETAG,QZQ,U12, &
            G12(3,3),A(3,3), &
            Zq(3), Z(3,3),Q12(3)
    real*8 :: UXY0(3,3), UX0(3), U

    if (y(3 * fm%Natom + 1) == 0d0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if

!    do I1=1,Natom-1
    call fm%get_g(fm%G) 

    U = 0; fm%UX = 0; fm%UXY = 0;
    do I1 = 1, fm%Natom - 1
        I1_3 = 3*(I1-1) + 1
        do I2 = I1+1, fm%Natom
            I2_3 = 3*(I2-1) + 1
            Q12 = y(I1_3:I1_3+2) - y(I2_3:I2_3+2)
            Q12 = min_image(Q12, fm%bl)
            G12=fm%G (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) &
                + fm%G (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) &
                - fm%G (I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) &
                - fm%G (I1_3 : I1_3 + 2, I2_3 : I2_3 + 2)

            call detminvm(G12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXY0 = 0d0
            DO IG=1, fm%NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
                AG = A
                do J=1,3
                    AG(J,J) = fm%LJA(IG) + AG(J,J)
                end do

                call detminvm(AG, DETAG, Z)
                Z = - fm%LJA(IG)**2 * Z

                do J=1,3
                    Z(J,J) = Z(J,J) + fm%LJA(IG)
                end do

                Zq = matmul(Z, Q12) ! R = -2.0*Zq
                qZq = dot_product(Q12, Zq) 

                U12 = SQRT(DETA/DETAG) * EXP(-qZq) * fm%LJC(IG)
                U = U + U12

                UX0 = UX0 - 2d0*U12*Zq
                do J=1,3
                    UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                end do
            end do ! IG
            
! Avoid load and store as much as possbile. Store now, process as a stream later. Much faster.
            fm%UX (I1_3 : I1_3 + 2) = fm%UX (I1_3 : I1_3 + 2) + UX0
            fm%UX (I2_3 : I2_3 + 2) = fm%UX (I2_3 : I2_3 + 2) - UX0

            fm%UXY (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) = &
                fm%UXY (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) + UXY0
            
            fm%UXY (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) = &
                fm%UXY (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) + UXY0
            
            !if (sum(q12**2) <= rfullmatsq) then
            fm%UXY (I1_3 : I1_3 + 2, I2_3 : I2_3 + 2) = -UXY0
            fm%UXY (I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) = -UXY0
            !end if
        end do ! I2
    end do ! I1


    !QP = -matmul(G, UX)
    call dsymm('L', 'L', 3 * fm%Natom, 1, -1d0, fm%G, 3 * fm%Natom, fm%UX, &
        3 * fm%Natom, 0d0, yp, 3 * fm%Natom)

    TRUXXG = sum( fm%UXY * fm%G)

    !GU = matmul(G, UXY)
    call dsymm('L', 'L', 3 * fm%Natom, 3 * fm%Natom, 1d0, fm%G, 3 * fm%Natom, &
        fm%UXY, 3 * fm%Natom, 0d0, fm%GU, 3 * fm%Natom)

    !GP = -matmul(GU, G)
    call dsymm('R', 'L', 3 * fm%Natom, 3 * fm%Natom, -1d0, fm%G, 3 * fm%Natom,&
        fm%GU, 3 * fm%Natom, 0d0, fm%GP, 3 * fm%Natom)

    do i1=1, fm%Natom
        k1=             3*(i1 - 1) + 1

        fm%GP (k1, k1)     = fm%GP (k1, k1)     + fm%invmass(i1)
        fm%GP (k1+1, k1+1) = fm%GP (k1+1, k1+1) + fm%invmass(i1)
        fm%GP (k1+2, k1+2) = fm%GP (k1+2, k1+2) + fm%invmass(i1)
    end do
    
    call fm_set_g(fm%GP, yp)

    yp(NEQ) = -(0.25d0*TRUXXG + U)/real(fm%Natom)
END SUBROUTINE RHSSFM

subroutine rhss_zero_time(NEQ, y, yp)
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), qi(3), qj(3), rsq, U
    integer :: i, j

    yp = 0d0

    j = 3 * fm%Natom + 1
    do i=1, 3 * fm%Natom
        yp(j) = fm%invmass
        j = j + 3 * fm%Natom - i + 1
    end do

    U=0d0

    DO I=1,fm%Natom-1
        qi = y(3*I-2:3*I)
        DO J=I+1,fm%Natom
                qj = y(3*J-2 : 3*J)
                qij = qi - qj
                rsq = sum(min_image(qij, fm%BL)**2)
                U = U + sum( fm%LJC (1 : fm%NGAUSS) &
                   * EXP( -fm%LJA (1 : fm%NGAUSS) * rsq) )
        ENDDO
    ENDDO

    yp(NEQ) = -U/real(fm%Natom)
end subroutine
end module vgwfm_mod
