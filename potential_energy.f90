subroutine Upot_tau0(Q,U)
    IMPLICIT NONE
    REAL*8, intent(in) :: Q(:,:)
    real*8, intent(out) :: U
    INTEGER  I,J,K, N
    real*8 :: RSQ,BL2,QIJ

    N = size(Q, 2)
    BL2=BL/2.0
    U=0.0

    DO I=1,N-1
        DO J=I+1,N
            RSQ=0.0D0
            DO K=1,3
                QIJ=Q(K,I)-Q(K,J)
                IF(QIJ > BL2) QIJ=QIJ-BL
                IF(QIJ < -BL2) QIJ=QIJ+BL
                RSQ=RSQ+QIJ**2
            ENDDO
            DO K=1,NGAUSS
                U=U+LJC(K)*EXP(-LJA(K)*RSQ)
            ENDDO
        ENDDO
    ENDDO
end subroutine Upot_tau0


function Ux_tau0(Q) result(UX)
    IMPLICIT NONE
    REAL*8, intent(in) :: Q(:,:)
    real*8 :: UX(3,size(Q,2))
    INTEGER  I,J,K,P, N
    REAL*8  RSQ,BL2,VEC(3),GRAD

    N = size(Q, 2)

    BL2=BL/2

    UX=0.0
    DO I=1,N-1
        DO J=I+1,N
            RSQ=0.0D0
            DO K=1,3
                VEC(K)=(Q(K,I)-Q(K,J))
                IF(VEC(K) > BL2) VEC(K)=VEC(K)-BL
                IF(VEC(K) < -BL2) VEC(K)=VEC(K)+BL
                RSQ=RSQ+VEC(K)**2
            ENDDO
            DO K=1,3
                GRAD=0.0D0
                DO P=1,NGAUSS
                    GRAD=GRAD+2.0*LJC(P)*LJA(P)*EXP(-RSQ*LJA(P))
                ENDDO
                GRAD=GRAD*VEC(K)
                UX(K,I)=UX(K,I)+GRAD
                UX(K,J)=UX(K,J)-GRAD
            ENDDO
        ENDDO
    ENDDO
end function Ux_tau0
