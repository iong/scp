SUBROUTINE RHSS1(NEQ, T, y, yprime)
      IMPLICIT NONE
    integer, intent(in) :: NEQ
      REAL*8, intent(in) :: T, y(NEQ)
      REAL*8, intent(out) :: yprime(NEQ)
      REAL*8 :: Q(3,N_atom), BLKC(3,3,N_atom), QNK(3,3,N_atom)
      INTEGER I,J,K,I1,I2,IG,CNT,CNT2
      REAL*8 AG(3,3),GU(3,3),&
           DETA,DETAG,GUG,QP,TRUXXGI,FACTOR,U,UX,UXX,QZQ,EXPAV, &
           TRMG,DETS,DETI,LJS,LJE,GUQ, &
           BL2,M(3,3),A(3,3), &
           R(3), Z(3,3),Q12(3)
      REAL*8 UPV(3,N_atom), UPM(3,3,N_atom)

    if (NEQ /= (1+21*N_atom) ) then
        write (*,*) 'NEQ != 1+21*N_atom', NEQ, 1+21*N_atom
        stop
    endif
    
    BL2=BL/2

    Q = reshape(y(2:1+3*N_atom), (/3, N_atom/) )
    call unpackg(N_atom, y(2+3*N_atom:1+9*N_atom), BLKC)

    CNT = 2+9*N_atom
    DO I=1,N_atom
        DO J=1,3
            DO K=1,3
                QNK(J,K,I)=Y(CNT)
                CNT=CNT+1
            ENDDO
        ENDDO
    ENDDO



    TRMG=0.0D0
    DO I=1,N_atom
        CALL INVDET(BLKC(:,:,I),M,DETS)
        DETS=1.0D0/DETS

        DO J=1,3
            TRMG=TRMG + M(J,J)*DETS
        ENDDO
    ENDDO

    TRMG = TRMG*MASS

    U=0.0D0
    UPV=0.0
    UPM = 0.0

      DO I1=1,N_atom-1
        DO I2=I1+1,N_atom
          IF(QRC(I2+(I1-1)*N_atom)) THEN
            DO I=1,3
              Q12(I)=Q(I,I1)-Q(I,I2)
              IF(Q12(I) > BL2) Q12(I)=Q12(I)-BL
              IF(Q12(I) < -BL2) Q12(I)=Q12(I)+BL
            ENDDO
            DO J=1,3
              DO K=J,3
                A(J,K)=BLKC(J,K,I1)+BLKC(J,K,I2)
              ENDDO
            ENDDO

            CALL INVDET(A,M,DETI)
            DETA=1.0D0/DETI

            DO J=1,3
              DO K=J,3
                A(J,K)=M(J,K)*DETA
              ENDDO
            ENDDO

            DO IG=1,NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
              DO J=1,3
                AG(J,J)=LJA(IG)+A(J,J)
                DO K=J+1,3
                  AG(J,K)=A(J,K)
                ENDDO
              ENDDO

              CALL INVDET(AG,M,DETAG)
              FACTOR=LJA(IG)**2/DETAG

              DO I=1,3
                Z(I,I)=LJA(IG)-FACTOR*M(I,I)
                DO J=I+1,3
                  Z(I,J)=-FACTOR*M(I,J)
                  Z(J,I)=Z(I,J)
                ENDDO
              ENDDO

              QZQ=0.0D0

              DO I=1,3
                R(I)=0.0D0
                DO J=1,3
                  R(I)=R(I)+Z(I,J)*Q12(J)
                ENDDO
                QZQ=QZQ+R(I)*Q12(I)
                R(I)=-2.0D0*R(I)
              ENDDO

              EXPAV=SQRT(DETA/DETAG)*EXP(-QZQ)
              U=U+EXPAV*LJC(IG)

              DO I=1,3
                UX=EXPAV*R(I)*LJC(IG)
                UPV(I,I1)=UPV(I,I1)+UX ! COMPUTE  GRADIENT (UPV = GRADIE
                UPV(I,I2)=UPV(I,I2)-UX
              DO J=I,3
                UXX=(R(I)*R(J)-2*Z(I,J))*EXPAV*LJC(IG) ! HESSIAN  MATRIX
                UPM(I,J,I1)=UPM(I,J,I1)+UXX ! UPM  = BLOCK  DIAGONAL OF
                UPM(I,J,I2)=UPM(I,J,I2)+UXX
                IF(I.NE.J) THEN ! FILL LOWER HALF OF MATRIX
                  UPM(J,I,I1)=UPM(J,I,I1)+UXX
                  UPM(J,I,I2)=UPM(J,I,I2)+UXX
                ENDIF
              ENDDO
            ENDDO
         ENDDO               ! END SUMMATION OVER GUASSIANS
       ENDIF
      ENDDO                  ! END OF I1 LOOP (ITH PARTICLE)
      ENDDO                     ! END OF I2 LOOP (JTH PARTICLE)

      TRUXXGI=0.0D0

      DO I1=1,N_atom
        DO I=1,3
          TRUXXGI=TRUXXGI+UPM(I,I,I1)*BLKC(I,I,I1)
          DO J=I+1,3
            TRUXXGI=TRUXXGI+2.0D0*UPM(I,J,I1)*BLKC(I,J,I1)
          ENDDO
        ENDDO
      ENDDO

       YPRIME(1)=-0.25D0*TRUXXGI-U
       CNT=2
       CNT2=2+9*N_atom

       DO I1=1,N_atom
         DO I=1,3
           QP=0.0D0
           DO J=1,3
             QP = QP-BLKC(I,J,I1)*UPV(J,I1)
           ENDDO
           YPRIME(CNT)=QP
           CNT=CNT+1
         ENDDO
       ENDDO
      
       DO I1=1,N_atom
         DO I=1,3
           DO J=1,3
             GU(I,J)=0.0D0
             DO K=1,3
               GU(I,J)=GU(I,J)+BLKC(I,K,I1)*UPM(K,J,I1)
             ENDDO
           ENDDO
        ENDDO
         
        DO I=1,3
          DO J=I,3
            GUG=0.0D0
            DO K=1,3
              GUG=GUG-GU(I,K)*BLKC(K,J,I1)
            ENDDO
            IF(I == J) THEN
              GUG=GUG+MASS
            ENDIF
            YPRIME(CNT)=GUG
            CNT=CNT+1
          ENDDO
        ENDDO

        DO I=1,3
          DO J=1,3
            GUQ=0.0D0
            DO K=1,3
              GUQ=GUQ-GU(I,K)*QNK(K,J,I1)
            ENDDO
            YPRIME(CNT2)=GUQ
            CNT2=CNT2+1
          ENDDO
        ENDDO
      ENDDO

      DO I=1,N_atom
        DO J=1,3
          QP=0.0D0
          DO K=1,3
            QP=QP-UPV(K,I)*QNK(K,J,I)
          ENDDO
          YPRIME(CNT2)=QP
          CNT2=CNT2+1
        ENDDO
      ENDDO
           

      return
END SUBROUTINE
