SUBROUTINE potential_energy(Q,N,U)
      use vgw
      IMPLICIT NONE
      INTEGER  I,J,K,N,P
      REAL*8  Q(3,N),U,RSQ,BL2,QIJ,VEC(3),GRAD
      
      BL2=BL/2
      U=0

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
END


