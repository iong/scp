      SUBROUTINE vgw0(N_ATOM,IMASS,QCNFG,W, TAUMAX,TAUI,BL,ATOL,RC,Y)
      IMPLICIT NONE
      integer LIW,MF,IERR
      real*8 RTOL,ATOL,ATOMICMASS,RC
      parameter(LIW=20,MF=10,ATOMICMASS=0.020614788876D0, &
           RTOL=0.0D0)
      INTEGER N_atom,NEQ,LRW,N_STEP,ITASK,ISTATE,IOPT,ITOL,I,J,l, &
           K,CNT,UF,IWORK(LIW)
      REAL*8 QCNFG(3,N_atom),ENRG,FX(3,N_atom),QFINAL(3,N_atom), &
           MASS,IMASS,W,LNZ(2),TAU(2),LOGZ,T,C0,UX(3,N_atom), &
           TOUT,ULJ,LNP,LJS,LJE,BL,BL2, &
           TAUMAX,TSTEP,TAUI,dummy,A(3,3),Z(3,3),F1(3),F2(3),eig(3), &
           Fscaled(3,N_atom),Fc(3)

      REAL*8 Y(1+9*N_atom), RWORK(36+336*N_atom)
      EXTERNAL RHSS0
      EXTERNAL JAC

      NEQ=1+9*N_atom
      LRW=36+336*N_atom
      TSTEP=0.1D0
      BL2=BL/2

      MASS=1.0D0/(ATOMICMASS*IMASS)

      call requal(3*N_atom,QCNFG,Y(2)) !Copy particle positions to Y

!      DO I=1,3*N_atom
!         Y(I+1)=Y(I+1)-BL*DNINT(Y(I+1)/BL) ! Periodic boundary conditions: wrap around box  ???????????????
!      ENDDO

      call initialize_SG_3G(MASS,BL,N_atom)
      call initialize_QRC(QCNFG,N_atom,RC,BL) !Determine which particles

      CALL DLSODEINIT(IWORK,RWORK,ITASK,IOPT,ISTATE,ITOL,LIW,LRW)

      C0=TAUI*MASS

      CNT=2+3*N_atom

      DO I=1,N_atom             ! INITIALIZE Y VECTOR WITH INITITIAL CON
        DO J=1,3
          Y(CNT)=C0
          CNT=CNT+1
          DO K=J+1,3
            Y(CNT)=0.0D0
            CNT=CNT+1
          ENDDO
        ENDDO
      ENDDO

      DO I=1,N_atom             ! INITIALIZE Q MATRIX
        DO J=1,3
          DO K=1,3
            IF(K.EQ.J) THEN
              Y(CNT)=1.0D0
            ELSE
              Y(CNT)=0.0D0
            ENDIF
            CNT=CNT+1
          ENDDO
        ENDDO
      ENDDO

      CALL GAUSSENERGYPB(QCNFG,N_atom,BL,ULJ,Y(CNT))
      CALL rscale(3*N_atom,Y(CNT),TAUI)
      Y(1)=-TAUI*ULJ

      T=TAUI
      IF(TAUMAX < 1) TSTEP = 0.001D0
      IF((TAUMAX/2)-TAUI < TSTEP) TSTEP = ((TAUMAX/2)-TAUI)/2.0D0

      TAU(1) = TAUMAX/2.0D0 - TSTEP
      TAU(2) = TAUMAX/2.0D0
!      write (*,*) y(1:4)
!      write(*,*) y(2+3*N_atom:7+3*N_atom),  y(2+9*N_atom:10+9*N_atom)
!      write (*,*) y(2+18*N_atom:4+18*N_atom)
        TOUT=TAU(2)
        CALL DLSODE(RHSS0,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE, &
                 IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
        CALL LNPS(N_atom,NEQ,Y,LOGZ)

      LNP=LOGZ
!        write (*,*) LOGZ, TAUMAX
      W=-(1/TAUMAX)*LNP-((1.5D0/TAUMAX)* &
              LOG(N_atom*2*3.1415926*TAUMAX/MASS))


      return
      END


      SUBROUTINE RHSS0(NEQ,T,Y,YPRIME)
      IMPLICIT NONE
      INTEGER NGAUSS,N_atom,I,J,K,I1,I2,IG,NEQ,CNT,CNT2,MAXN
      PARAMETER (NGAUSS=4,MAXN=5000)
      REAL*8 T,AG(3,3),GU(3,3),Y(NEQ+3*(NEQ-1)/9),YPRIME(NEQ), &
           DETA,DETAG,GUG,QP,TRUXXGI,FACTOR,U,UX,UXX,QZQ,EXPAV, &
           TRMG,DETS,DETI,MASS,LJS,LJE,GUQ, &
           BL,BL2,M(3,3),A(3,3), &
           R(3), Z(3,3),Q12(3)
      REAL*8 Q(3,(NEQ-1)/9),UPV(3,(NEQ-1)/9),BLKC(3,3,(NEQ-1)/9), &
           UPM(3,3,(NEQ-1)/9)

      real*8 LJC(NGAUSS),LJA(NGAUSS)
      logical QRC(MAXN*MAXN)
      common /LJ3/ LJC,LJA
      common /systemee/ MASS,BL,N_atom
      common /carray/ QRC
!     equivalence (Y(1),gamma), (Y(2),Q(1,1))

!      write (*,*) T
      if(N_atom.ne.(NEQ-1)/9) stop 'RHSS0: Oops!'
      BL2=BL/2

      call requal(3*N_atom,Y(2),Q)
      TRMG=0.0D0

      CNT=3*N_atom+2

      DO I=1,N_atom
        DO J=1,3
          DO K=J,3
            BLKC(J,K,I)=Y(CNT)
            BLKC(K,J,I)=Y(CNT)
            CNT=CNT+1
          ENDDO
        ENDDO
      ENDDO

      DO I=1,N_atom
        call requal(9,BLKC(1,1,I),A)
        CALL INVDET(A,M,DETS)
        DETS=1.0D0/DETS

        DO J=1,3
          TRMG=TRMG + M(J,J)*DETS
        ENDDO
      ENDDO

      TRMG = TRMG*MASS

      U=0.0D0
      call rzero(3*N_atom,UPV)
      call rzero(9*N_atom,UPM)

      DO I1=1,N_atom-1
        DO I2=I1+1,N_atom
          IF(QRC(I2+(I1-1)*N_atom)) THEN
            DO I=1,3
              Q12(I)=Q(I,I1)-Q(I,I2)
              IF(Q12(I) > BL2) Q12(I)=Q12(I)-BL      !  fix this later
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

       CNT=1
       YPRIME(CNT)=-0.25D0*TRUXXGI-U
       CNT=CNT+1
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
      ENDDO

      return
      END
