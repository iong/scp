SUBROUTINE Upot_tau0(Q,N,U)
    use vgw
    IMPLICIT NONE
    integer, intent(in) :: N
    REAL*8, intent(in) :: Q(3,N)
    real*8, intent(out) :: U
    INTEGER  I,J,K,P
    real*8 :: RSQ,BL2,QIJ,VEC(3),GRAD
      
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
      return
END SUBROUTINE


SUBROUTINE Ux_tau0(Q,N,UX)
    use vgw
    IMPLICIT NONE
    integer, intent(in) :: N
    REAL*8, intent(in) :: Q(3,N)
    real*8, intent(out) :: UX(3,N)
    INTEGER  I,J,K,P
    REAL*8  RSQ,BL2,QIJ,VEC(3),GRAD
      
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
      return
END SUBROUTINE
