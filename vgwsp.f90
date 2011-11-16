module vgwsp_mod
    use vgw_mod
    implicit none

    type, public, extends(vgw) :: vgwsp
        real*8, allocatable :: UPV(:,:), UPM(:,:,:), G(:,:,:), Qnk(:,:,:)
    contains
        procedure :: init => vgwsp_init
        procedure :: cleanup => vgwsp_cleanup
        procedure :: init_Qnk0 => vgwsp_init_Qnk0
        procedure :: propagate => vgwsp_propagate
        procedure :: logdet => vgwsp_logdet
    end type vgwsp

    type(vgwsp) :: sp

contains

subroutine vgwsp_init(this, Natom, species)
    implicit none
    class(vgwsp) :: this
    integer, intent(in) :: Natom
    character(*), intent(in) :: species

    this%Natom = Natom
    this%NG = 6*Natom
    this%NQnk = 9*Natom

    this%species = species

    allocate(this%G(3,3,Natom), this%Qnk(3,3,Natom), this%UPV(3,Natom), &
        this%UPM(3,3,Natom) )
    
    call this%common_init()
end subroutine


subroutine vgwsp_cleanup(this)
    implicit none
    class(vgwsp) :: this

    deallocate(this%upv, this%upm, this%g, this%Qnk)

    call this%common_destroy()
end subroutine


subroutine sp_get_g(N, y, g)
    implicit none
    integer, intent(in) :: N
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: g(:,:,:)
    integer :: i, k

    i=0
    do k=3*N+1,9*N,6
        i=i+1
        g(:,1,i) = y(k : k+2)
        g(1, 2, i) = y(k+1)
        g(2:3,2,i) = y(k+3 : k+4)
        g(1, 3, i) = y(k+2)
        g(2, 3, i) = y(k+4)
        g(3, 3, i) = y(k+5)
    end do
end subroutine


subroutine sp_get_Qnk(N, y, Qnk)
    implicit none
    integer, intent(in) :: N
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: Qnk(:,:,:)

    Qnk = reshape(y(9*N+1 : 18*N), (/3, 3, N/))
end subroutine


subroutine vgwsp_init_Qnk0(this)
    implicit none
    class(vgwsp) :: this
    integer :: j

    do j=1, this%Natom
        this%y(9*this%Natom + 9*(j - 1) + 1) = 1d0
        this%y(9*this%Natom + 9*(j - 1) + 5) = 1d0
        this%y(9*this%Natom + 9* j         ) = 1d0
    end do
end subroutine


function vgwsp_logdet(this)
    implicit none
    class(vgwsp) :: this
    double precision :: vgwsp_logdet

    integer :: j, Gi

    vgwsp_LOGDET=0d0
    DO j=1,this%Natom
        Gi = 3 * this%Natom + 6*(j-1) + 1
        vgwsp_LOGDET = vgwsp_LOGDET + LOG( DETM_S( this%y(Gi : Gi+5) ) )
    ENDDO
end function

subroutine vgwsp_propagate(this, tau)
    implicit none
    class(vgwsp) :: this
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
            DETA,DETAG,GUG(6),UX(3,sp%nnbmax),UXX(3,3, sp%nnbmax),QZQ,EXPAV, &
            G12(6),A(3,3), &
            Zq(3), Z(3,3),Q12(3), v0, QP_(3), Q1(3), G1(6)
    real*8 :: GC(6,sp%nnbmax), &
            QC(3, sp%nnbmax), UXX0(3,3), UX0(3), UPV_I1(3), UPM_I1(3,3)
    real*8 :: U, TRUXXG

    sp%UPM = 0d0
    sp%UPV = 0d0
    U = 0d0

    if (y(NEQ) == 0d0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if
        

    do I1 = 1, sp%Natom - 1
        NN1 = sp%NNB(I1)
        if (NN1 == 0) cycle

        q1 = y(3*I1 - 2 : 3*I1)
        G1 = y(3 * sp%Natom + 6*I1 - 5 : 3 * sp%Natom + 6*I1)
        UPV_I1 = sp%UPV (:,I1) 
        UPM_I1 = sp%UPM (:,:,I1) 
        
        do I2 = 1, NN1
            k1 = 3 * sp%nbidx(I2,I1) - 2
            k2 = 3 * sp%Natom + 6 * sp%nbidx(I2,I1) - 5
            qc(:,I2) = y(k1 : k1 + 2)
            GC(:,I2) = y(k2 : k2 + 5)
        end do

        DO I2=1,NN1
            Q12 = Q1 - QC(:,I2)
            Q12 = min_image(Q12, sp%bl)
            G12=G1+GC(:,I2)

            call detminvm_sg(G12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXX0 = 0d0
            DO IG=1, sp%NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
                AG = A
                do J=1,3
                    AG(J,J)=sp%LJA(IG)+AG(J,J)
                end do

                call detminvm(AG, DETAG, Z)
                Z = - sp%LJA(IG)**2 * Z

                do J=1,3
                    Z(J,J) = Z(J,J) + sp%LJA(IG)
                end do

                Zq = matmul(Z, Q12) ! R = -2.0*Zq
                qZq = dot_product(Q12, Zq) 

                EXPAV=SQRT(DETA/DETAG)*EXP(-qZq)
                U = U + EXPAV * sp%LJC(IG)

                v0 = 2d0 * expav * sp%LJC(IG)

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
        sp%UPV (:,I1) = UPV_I1

        sp%UPM (:,:,I1) = UPM_I1
        
        sp%UPV (:,sp%nbidx(1:NN1,I1)) = sp%UPV (:,sp%nbidx(1:NN1,I1)) &
            + UX(:,1:NN1)
        
        sp%UPM (:,:,sp%nbidx(1:NN1,I1)) = sp%UPM (:,:,sp%nbidx(1:NN1,I1)) &
            + UXX(:,:,1:NN1)
    ENDDO ! I

    call sp_get_g(sp%Natom, sp%y, sp%g)

    if (NEQ > 21*sp%Natom ) then
        call sp_get_Qnk (sp%Natom, sp%y, sp%Qnk)
    end if

    TRUXXG = sum(sp%UPM * sp%G)

    do i1=1,sp%Natom
        k1=             3*(i1 - 1) + 1
        k2 =  3*sp%Natom + 6*(i1 - 1) + 1
        k3 =  9*sp%Natom + 9*(i1 - 1) + 1
        k4 = 18*sp%Natom + 3*(i1 - 1) + 1

        QP_ = - matmul(sp%G(:,:,I1), sp%UPV (:,I1))

        GU = matmul(sp%G(:,:,I1), sp%UPM (:,:,I1))
        GUG = -matmul_sgs(GU, sp%G (:,:,I1) )
        GUG(1) = GUG(1) + sp%invmass (i1)
        GUG(4) = GUG(4) + sp%invmass (i1)
        GUG(6) = GUG(6) + sp%invmass (i1)

        yp(k1 : k1+2) = QP_
        yp(k2 : k2+5) = GUG

        if (NEQ > 21*sp%Natom ) then
            yp(k3 : k3+8) = - reshape(matmul(GU, sp%Qnk (:,:,i1)), (/ 9 /) )

            yp(k4 : k4 + 2 ) = - matmul(transpose(sp%Qnk (:,:,i1)), sp%UPV (:,I1))
        end if
    end do

    yp(NEQ) = -(0.25d0*TRUXXG + U)/real(sp%Natom)
END SUBROUTINE

function matmul_sgs(A, B) result (C)
    implicit none
    double precision :: A(3,3), B(3,3)
    double precision :: C(6)

    C(1:3) = matmul(A, B(:,1))
    C(4) = sum(A(2,:)*B(:,2))
    C(5) = sum(A(3,:)*B(:,2))
    C(6) = sum(A(3,:)*B(:,3))
end function

subroutine rhss_zero_time(NEQ, y, yp)
    implicit none
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), qi(3), qj(3), rsq, LJC_EXPAV(sp%NGAUSS), &
        UX0(3), U
    integer :: i, j, ki, kj

    yp = 0d0

    do i=1,sp%Natom
        ki = 3*sp%Natom + 6*i - 5
        yp(ki : ki + 5) = sp%invmass(i) * (/1d0, 0d0, 0d0, 1d0, 0d0, 1d0/)
    end do

    U=0d0

    
    DO I=1, sp%Natom-1
        qi = y(3*I-2:3*I)
        DO J=1, sp%NNB(I)
                qj = y(3 * sp%NBIDX(J,I)-2 : 3 * sp%NBIDX(J,I))
                qij = qi - qj
                rsq = sum(min_image(qij, sp%BL)**2)

                LJC_EXPAV = sp%LJC (1 : sp%NGAUSS) * EXP(-rsq &
                    * sp%LJA (1 : sp%NGAUSS))

                U = U + sum(LJC_EXPAV)

                ! there are to -1 here which cancel each other
                UX0 = qij*(2.0*sum( sp%LJA (1 : sp%NGAUSS) * LJC_EXPAV ) )

                ki = 3*sp%Natom + sp%NG + sp%NQnk + 3*I  - 2
                kj = 3*sp%Natom + sp%NG + sp%NQnk + 3*J  - 2

                yp(ki : ki + 2) = yp(ki : ki + 2) + UX0
                yp(kj : kj + 2) = yp(kj : kj + 2) - UX0
        ENDDO
    ENDDO

    yp(NEQ) = -U/real(sp%Natom)
end subroutine

end module vgwsp_mod
