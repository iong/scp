module vgwsp_mod
    use vgw_mod
    implicit none

    type, public, extends(vgw) :: vgwsp
        real*8, allocatable :: UPV(:,:), UPM(:,:,:), G(:,:,:), Qnk(:,:,:)
    contains
        procedure :: init => vgwsp_init
        procedure :: Ueff => vgwsp_Ueff
        procedure :: cleanup => vgwsp_cleanup
    end type vgwsp

contains

subroutine vgwsp_init(Nmax_, species, rcutoff, massx)
    implicit none
    integer, intent(in) :: Nmax_
    character(*), intent(in) :: species
    real*8, intent(in), optional :: rcutoff, massx(:)

    Nmax = Nmax_
    NG = 6*Nmax
    NQNK = 9*Nmax
    NEQMAX = 3*Nmax + NG + NQNK + 3*Nmax  + 1

    allocate(NNB(Nmax), NBIDX(Nmax,Nmax), upv(3,Nmax), upm(3,3,Nmax), &
        g(3,3,Nmax), Qnk(3,3,Nmax), mass(Nmax), invmass(Nmax), y(NEQMAX))
    
include 'species.f90'
end subroutine


subroutine vgwsp_cleanup()
    deallocate(NNB, NBIDX, UPV, UPM, g, Qnk, mass, invmass, y)
end subroutine


subroutine sp_get_g(y, g)
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: g(:,:,:)
    integer :: i, k

    i=0
    do k=3*Natom+1,9*Natom,6
        i=i+1
        g(:,1,i) = y(k : k+2)
        g(1, 2, i) = y(k+1)
        g(2:3,2,i) = y(k+3 : k+4)
        g(1, 3, i) = y(k+2)
        g(2, 3, i) = y(k+4)
        g(3, 3, i) = y(k+5)
    end do
end subroutine


subroutine sp_get_Qnk(y, qnk)
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: qnk(:,:,:)

    q = reshape(y(9*Natom+1 : 18*Natom), (/3, 3, Natom/))
end subroutine


subroutine sp_init_Qnk0(this)
    class(vgwfm) :: this
    integer :: j

    do j=1,Natom
        this%y(9*this%Natom + 9*(j - 1) + 1) = 1d0
        this%y(9*this%Natom + 9*(j - 1) + 5) = 1d0
        this%y(9*this%Natom + 9* j         ) = 1d0
    end do
end subroutine


function sp_logdet(this)
    class(vgwsp) :: this

    sp_LOGDET=0d0
    DO j=1,Natom
        Gi = 3 * this%Natom + 6*(j-1) + 1
        sp_LOGDET = sp_LOGDET + LOG( DETM_S( this%y(Gi : GI+5) ) )
    ENDDO
end function sp_logdet

subroutine vgwsp_propagate(this, tau)
    class(vgwfm) :: this
    double precision :: tau

    call this%dlsode%propagate(RHSSP, this%y, this%T, tau)
end subroutine



SUBROUTINE RHSSP(NEQ, T, Y, YP)
    IMPLICIT NONE
    integer, intent(in) :: NEQ!, IPAR(:)
    double precision, intent(in) :: T, Y(NEQ)!, RPAR(:)
    double precision, intent(out) :: YP(NEQ)

    INTEGER :: J,I1,I2,IG, NN1, k1, k2, k3, k4
    REAL*8 AG(3,3),GU(3,3), &
            DETA,DETAG,GUG(6),UX(3,nnbmax),UXX(3,3, nnbmax),QZQ,EXPAV, &
            G12(6),A(3,3), &
            Zq(3), Z(3,3),Q12(3), v0, QP_(3), Q1(3), G1(6)
    real*8 :: GC(6,nnbmax), &
            QC(3, nnbmax), UXX0(3,3), UX0(3), UPV_I1(3), UPM_I1(3,3)
    real*8 :: Ulocal

    UPM = 0d0
    UPV = 0d0
    Ulocal = 0d0

    if (y(3*Natom+1) == 0d0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if
        

    do I1=1,Natom-1
        NN1 = NNB(I1)
        if (NN1 == 0) cycle

        q1 = y(3*I1 - 2 : 3*I1)
        G1 = y(3*Natom + 6*I1 - 5 : 3*Natom + 6*I1)
        UPV_I1 = UPV(:,I1) 
        UPM_I1 = UPM(:,:,I1) 
        
        do I2=1,NN1
            qc(:,I2) = y(3*nbidx(I2,I1) - 2 : 3*nbidx(I2,I1))
            GC(:,I2) = y(3*Natom + 6*nbidx(I2,I1) - 5 : 3*Natom + 6*nbidx(I2,I1))
        end do

        DO I2=1,NN1
            Q12 = Q1 - QC(:,I2)
            Q12 = min_image(Q12, bl)
            G12=G1+GC(:,I2)

            call detminvm_sg(G12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXX0 = 0d0
            DO IG=1,NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
                AG = A
                do J=1,3
                    AG(J,J)=LJA(IG)+AG(J,J)
                end do

                call detminvm(AG, DETAG, Z)
                Z = - LJA(IG)**2 * Z

                do J=1,3
                    Z(J,J) = Z(J,J) + LJA(IG)
                end do

                Zq = matmul(Z, Q12) ! R = -2.0*Zq
                qZq = dot_product(Q12, Zq) 

                EXPAV=SQRT(DETA/DETAG)*EXP(-qZq)
                Ulocal=Ulocal+EXPAV*LJC(IG)

                v0 = 2d0*expav*LJC(IG)

                UX0 = UX0 - v0*Zq
                do J=1,3
                    UXX0(:,J) = UXX0(:,J) + v0*(2d0*Zq*Zq(J) - Z(:,J))
                end do
            ENDDO ! IG
! Avoid load and store as much as possbile. Store now, process as a stream later. Much faster.
            UPV_I1 = UPV_I1 + UX0
            UX(:,I2) = - UX0
            UPM_I1 = UPM_I1 + UXX0
            UXX(:,:,I2) = UXX0
        ENDDO ! I2
        UPV(:,I1) = UPV_I1
        UPM(:,:,I1) = UPM_I1
        UPV(:,nbidx(1:NN1,I1)) = UPV(:,nbidx(1:NN1,I1)) + UX(:,1:NN1)
        UPM(:,:,nbidx(1:NN1,I1)) = UPM(:,:,nbidx(1:NN1,I1)) + UXX(:,:,1:NN1)
    ENDDO ! I

    call get_g(y, g)

    if (NEQ > 21*Natom ) then
        call get_Qnk(y, Qnk)
    end if

    TRUXXG = sum(UPM*G)
    U = Ulocal

    do i1=1,Natom
        k1=             3*(i1 - 1) + 1
        k2 =  3*Natom + 6*(i1 - 1) + 1
        k3 =  9*Natom + 9*(i1 - 1) + 1
        k4 = 18*Natom + 3*(i1 - 1) + 1

        QP_ = - matmul(G(:,:,I1), UPV(:,I1))

        GU = matmul(G(:,:,I1), UPM(:,:,I1))
        GUG = -matmul_sgs(GU, G(:,:,I1))
        GUG(1) = GUG(1) + invmass(i1)
        GUG(4) = GUG(4) + invmass(i1)
        GUG(6) = GUG(6) + invmass(i1)

        yp(k1 : k1+2) = QP_
        yp(k2 : k2+5) = GUG

        if (NEQ > 21*Natom ) then
            yp(k3 : k3+8) = - reshape(matmul(GU, Qnk(:,:,i1)), (/ 9 /) )

            yp(k4 : k4 + 2 ) = - matmul(transpose(Qnk(:,:,i1)), UPV(:,I1))
        end if
    end do

    yp(NEQ) = -(0.25d0*TRUXXG + U)/real(Natom)
END SUBROUTINE RHSS0

function matmul_sgs(A, B) result (C)
    double precision :: A(3,3), B(3,3)
    double precision :: C(6)

    C(1:3) = matmul(A, B(:,1))
    C(4) = sum(A(2,:)*B(:,2))
    C(5) = sum(A(3,:)*B(:,2))
    C(6) = sum(A(3,:)*B(:,3))
end function

subroutine rhss_zero_time(NEQ, y, yp)
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), qi(3), qj(3), rsq, LJC_EXPAV(NGAUSS), UX0(3)
    integer :: i, j

    yp = 0d0

    do i=1,Natom
        yp(3*Natom + 6*i - 5 : 3*Natom + 6*i) = invmass(i)*(/1d0, 0d0, 0d0, 1d0, 0d0, 1d0/)
    end do

    U=0d0

    DO I=1,Natom-1
        qi = y(3*I-2:3*I)
        DO J=1,NNB(I)
                qj = y(3*NBIDX(J,I)-2 : 3*NBIDX(J,I))
                qij = qi - qj
                rsq = sum(min_image(qij, BL)**2)
                LJC_EXPAV=LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq)
                U = U + sum(LJC_EXPAV)
                UX0 = (-1.0) * (-1.0) * qij*(2.0*sum(LJA(1:NGAUSS)*LJC_EXPAV))

                yp(18*Natom + 3*(I-1) : 18*Natom + 3*I) = &
                    yp(18*Natom + 3*(I-1) : 18*Natom + 3*I) + UX0

                yp(18*Natom + 3*(J-1) : 18*Natom + 3*J) = &
                    yp(18*Natom + 3*(J-1) : 18*Natom + 3*J) - UX0
                
        ENDDO
    ENDDO

    yp(NEQ) = -U/real(Natom)
end subroutine

