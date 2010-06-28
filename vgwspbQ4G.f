      SUBROUTINE vgwquenchspb(N_ATOM,IMASS,QCNFG,FX,LNP,W,
     &     ENRG,TAUMAX,TAUI,BL,RETQ)
      IMPLICIT NONE
      integer LIW,MF,RETQ
      real*8 RTOL,ATOL,ATOMICMASS
      parameter(LIW=20,MF=10,ATOMICMASS=0.020614788876D0,
     &     RTOL=0.0D0,ATOL=1.0D-3)
      INTEGER N_atom,NEQ,LRW,N_STEP,ITASK,ISTATE,IOPT,ITOL,I,J,
     &     K,CNT,UF,IWORK(LIW)
      REAL*8 QCNFG(3,N_atom),ENRG,FX(3,N_atom),QFINAL(3,N_atom),
     &     MASS,IMASS,W,LNZ(2),TAU(2),LOGZ,T,C0,UX(3,N_atom),
     &     TOUT,ULJ,LNP,LJS,LJE,BL,BL2,
     &     TAUMAX,TSTEP,TAUI

      REAL*8 Y(1+21*N_atom), RWORK(36+336*N_atom)

      EXTERNAL RHSS
      EXTERNAL RHSSIG
      EXTERNAL JAC
      
      NEQ=1+21*N_atom
      LRW=36+336*N_atom
      TSTEP=0.1D0
      BL2=BL/2
          
      MASS=1.0D0/(ATOMICMASS*IMASS)

      call initialize_SG_3G(MASS,BL,N_atom)
         
      CALL DLSODEINIT(IWORK,RWORK,ITASK,IOPT,ISTATE,ITOL,LIW,LRW)

      C0=TAUI*MASS
      call requal(3*N_atom,QCNFG,Y(2))
      
      CNT=2+3*N_atom
      
      DO I=1,N_atom             ! INITIALIZE Y VECTOR WITH INITITIAL CONDITIONS
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
     
      DO I=1,3*N_atom
        Y(I+1)=Y(I+1)-BL*INT(Y(I+1)/BL2)    ! Periodic boundary conditions: wrap around box
      ENDDO

      TAU(1) = TAUMAX/2.0D0 - TSTEP
      TAU(2) = TAUMAX/2.0D0
     
      DO I=1,2
        TOUT=TAU(I)
        CALL DLSODE(RHSS,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     &           IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
        CALL LNPS(N_atom,NEQ,Y,LOGZ)
        LNZ(I) = LOGZ+1.5D0*LOG(TOUT)
      ENDDO
  
      ENRG=-(LNZ(2)-LNZ(1))/(2*TSTEP)
      LNP=LOGZ
      W=-(1/TAUMAX)*LNP-((1.5D0/TAUMAX)*
     &        LOG(N_atom*2*3.1415926*TAUMAX/MASS))
    
      IF(RETQ==1) call requal(3*N_atom,Y(2),QFINAL)

      call requal(3*N_atom,Y(2+18*N_atom),FX)
      call rscale(3*N_atom,FX,1/TAU(2))
              
      IF(RETQ==1) call requal(3*N_atom,QFINAL,QCNFG)
  
      return
      END 
      
      SUBROUTINE DLSODEINIT(IWORK,RWORK,ITASK,IOPT,ISTATE,ITOL,LIW,LRW)
      IMPLICIT NONE
      INTEGER  ITASK,IOPT,ISTATE,LIW,LRW,ITOL,IWORK(LIW)
      REAL*8  RWORK(LRW)
      
      ITOL=1
      ITASK=1
      IOPT=1
      ISTATE=1
      IWORK(6)=100000
      IWORK(5)=4
      IWORK(7)=0
      IWORK(8)=0
      IWORK(9)=0
      IWORK(10)=0
      RWORK(5)=0.0D0
      RWORK(6)=0.0D0
      RWORK(7)=0.0D0
      RWORK(8)=0.0D0
      RWORK(9)=0.0D0
      RWORK(10)=0.0D0
      return
      END
      
      SUBROUTINE JAC 
      return
      END
      
      SUBROUTINE LNPS(N_atom,NEQ,Y,LOGZ)
      IMPLICIT NONE
      INTEGER  I,J,K,N_atom,NEQ,CNT
      real*8  Y(NEQ), C(3,3), M(3,3), LOGZ, GAMMA, DETI, DET
      
      EXTERNAL INVDET
      
      GAMMA=Y(1)
      CNT=2+3*N_atom
      
      DET=0.0D0
      
      DO I=1,N_atom
        DO J=1,3
          C(J,J) = Y(CNT)
          CNT=CNT+1
          DO K=J+1,3
            C(J,K)= Y(CNT)
            C(K,J)= C(J,K)
            CNT=CNT+1
          ENDDO
        ENDDO
        CALL INVDET(C,M,DETI)
        DET = DET + LOG(DETI)
      ENDDO
      
      LOGZ = 2.0D0*GAMMA - 0.5D0*DET
      return
      END
      
      SUBROUTINE INVDET(A, M, DET)
      IMPLICIT NONE
      REAL*8  DET, A(3,3), M(3,3)
      
      M(1,1) = A(2,2)*A(3,3)-A(2,3)*A(2,3)
      M(2,2) = A(1,1)*A(3,3)-A(1,3)*A(1,3)
      M(3,3) = A(1,1)*A(2,2)-A(1,2)*A(1,2)
      M(1,2) = -A(1,2)*A(3,3)+A(1,3)*A(2,3)
      M(1,3) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
      M(2,3) = -A(1,1)*A(2,3)+A(1,3)*A(1,2)
      DET = M(1,1)*A(1,1)+M(1,2)*A(1,2)+M(1,3)*A(1,3)
      return
      END 
      
      SUBROUTINE RHSS(NEQ,T,Y,YPRIME)
      IMPLICIT NONE
      EXTERNAL RHSSIG
      INTEGER NGAUSS,N_atom,I,J,K,I1,I2,IG,NEQ,CNT,CNT2
      PARAMETER  (NGAUSS=4)
      REAL*8 T,AG(3,3),GU(3,3),Y(NEQ+3*(NEQ-1)/21),YPRIME(NEQ),
     &     DETA,DETAG,GUG,QP,TRUXXGI,FACTOR,U,UX,UXX,QZQ,EXPAV,
     &     TRMG,DETS,DETI,MASS,LJS,LJE,GUQ,
     &     BL,BL2,M(3,3),A(3,3), 
     &     R(3), Z(3,3),Q12(3)
      REAL*8 Q(3,(NEQ-1)/21),UPV(3,(NEQ-1)/21),BLKC(3,3,(NEQ-1)/21),
     &     UPM(3,3,(NEQ-1)/21), QNK(3,3,(NEQ-1)/21) 
      REAL*8 YP2(NEQ)
      
      real*8 LJC(NGAUSS),LJA(NGAUSS)
      common /LJ3/ LJC,LJA
      common /systemee/ MASS,BL,N_atom
c     equivalence (Y(1),gamma), (Y(2),Q(1,1))

      if(N_atom.ne.(NEQ-1)/21) stop 'RHSS: Oops!'
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
        DO J=1,3
          DO K=1,3
            QNK(J,K,I)=Y(CNT)
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
              UPV(I,I1)=UPV(I,I1)+UX ! COMPUTE  GRADIENT (UPV = GRADIENT VECTOR)
              UPV(I,I2)=UPV(I,I2)-UX
            DO J=I,3
              UXX=(R(I)*R(J)-2*Z(I,J))*EXPAV*LJC(IG) ! HESSIAN  MATRIX (BLOCK OF I,JTH PARTICLE)
              UPM(I,J,I1)=UPM(I,J,I1)+UXX ! UPM  = BLOCK  DIAGONAL OF 3Nx3N HESSIAN
              UPM(I,J,I2)=UPM(I,J,I2)+UXX
              IF(I.NE.J) THEN ! FILL LOWER HALF OF MATRIX
                UPM(J,I,I1)=UPM(J,I,I1)+UXX
                UPM(J,I,I2)=UPM(J,I,I2)+UXX
              ENDIF
            ENDDO
          ENDDO
       ENDDO               ! END SUMMATION OVER GUASSIANS
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

c      CALL RHSSIG(NEQ,  T, Y, YP2, YPRIME)

      return
      END
      
      SUBROUTINE CLRADIUS(Q, N, RAD)
      IMPLICIT NONE
      INTEGER  N, I, J
      REAL*8  Q(3,N), RCM(3), CV(3), CR, RAD
      
      RAD=0
      call rzero(3,RCM)
      
      DO I=1,N
         DO J=1,3
            RCM(J)= RCM(J)+Q(J,I)
         ENDDO
      ENDDO
      
      DO I=1,3
         RCM(I)=RCM(I)/N
      ENDDO
      
      DO I=1,N
         CR=0
         DO J=1,3
            CV(J)=Q(J,I)-RCM(J)
            CR=CR+CV(J)**2
         ENDDO
         CR = SQRT(CR)
         IF(CR > RAD) THEN
            RAD=CR
         ENDIF
      ENDDO
      return
      END

      subroutine initialize_SG_3G(MASS,BL,N_atom)
      implicit none
      integer NGAUSS,N_atom,N_atom1
      parameter(NGAUSS=4)
      real*8 LJC(NGAUSS),LJA(NGAUSS),MASS,BL,
     &       MASS1,BL1
      common /LJ3/ LJC,LJA
      common /systemee/ MASS1,BL1,N_atom1
      data LJC /96609.488289873d0, 14584.62075507514d0, 
     &     -365.460614956589d0, -19.5534697800036d0/
     &     LJA /1.038252215127D0, 0.5974039109464D0,
     &     0.196476572277834D0, 0.06668611771781D0/
      logical initializeee
      data initializeee/.true./
      if (initializeee) then
        initializeee = .false.
      else
         return
      endif
      write(*,*) 'Initializing the 3 Gaussian LJ-potential'
      MASS1=MASS
      BL1=BL
      N_atom1=N_atom
      return
      end

      SUBROUTINE GAUSSENERGYPB(Q,N,BL,U,UX)
      IMPLICIT NONE
      INTEGER  I,J,K,N,P
      REAL*8  Q(3,N),U,RSQ,BL,BL2,QIJ,UX(3,N),VEC(3),GRAD
      integer NGAUSS
      parameter(NGAUSS=4)
      real*8 LJC(NGAUSS),LJA(NGAUSS)
      common /LJ3/ LJC,LJA
      
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

      call rzero(3*N,UX)
          
      DO I=1,N-1
        DO J=I+1,N
          RSQ=0.0D0
          DO K=1,3
            VEC(K)=(Q(K,I)-Q(K,J))
            IF(VEC(K).GT.BL2) VEC(K)=VEC(K)-BL
            IF(VEC(K).LT.-BL2) VEC(K)=VEC(K)+BL
            RSQ=RSQ+VEC(K)**2
          ENDDO
          DO K=1,3
            GRAD=0.0D0
            DO P=1,NGAUSS
              GRAD=GRAD+2*LJC(P)*LJA(P)*EXP(-RSQ*LJA(P))
            ENDDO
            GRAD=GRAD*VEC(K)
            UX(K,I)=UX(K,I)+GRAD
            UX(K,J)=UX(K,J)-GRAD
          ENDDO
        ENDDO
      ENDDO
      return
      END
      
      SUBROUTINE FORCES(N_atom,NEQ,Y,FX,Q,TAU)
      IMPLICIT NONE
      INTEGER I,J,K,IND,N_atom,NEQ
      REAL*8  TAU,Y(NEQ),FX(3,N_atom),Q(3,N_atom),QX(3,N_atom),
     &     BLK(3,3),INV(3,3)
      
      EXTERNAL INVS
      
      IND=2
      call rzero(9,BLK)
      call rzero(3*N_atom,FX)
      
      DO I=1,N_atom
         do J=1,3
            QX(J,I)=Y(IND)-Q(J,I)
            IND=IND+1
         enddo
      ENDDO
      
      DO I=1,N_atom
         DO J=1,3
            DO K=J,3
               BLK(K,J)=Y(IND)
               BLK(J,K)=Y(IND)
               IND=IND+1
            ENDDO
         ENDDO
         CALL INVS(BLK,INV)
         DO J=1,3
            DO K=1,3
               FX(J,I)=FX(J,I)+INV(K,J)*QX(K,I)
            ENDDO
            FX(J,I)=FX(J,I)*(2/TAU)
         ENDDO
      ENDDO
      return
      END
      subroutine rzero(N,A)
      implicit none
      integer i,N
      real*8 A(N)
      do i=1,N
         A(i)=0
      enddo
      return
      end

      subroutine requal(N,A,B)
      implicit none
      integer i,N
      real*8 A(N),B(N)
      do i=1,N
         B(i)=A(i)
      enddo
      return
      end

      subroutine rscale(N,A,c)
      implicit none
      integer i,N
      real*8 A(N),c
      do i=1,N
         A(i)=A(i)*c
      enddo
      return
      end
          
      SUBROUTINE INVS(M,MI)
      IMPLICIT NONE
      REAL*8  DET,M(3,3),MI(3,3)
      EXTERNAL INVDET
      
      DET=M(1,1)*(M(2,2)*M(3,3)-M(3,2)*M(2,3))-M(1,2)*(M(2,1)*M(3,3)-
     &     M(3,1)*M(2,3))+M(1,3)*(M(1,2)*M(3,2)-M(3,1)*M(2,2))
      
      MI(1,1) = M(2,2)*M(3,3)-M(2,3)*M(2,3)
      MI(2,2) = M(1,1)*M(3,3)-M(1,3)*M(1,3)
      MI(3,3) = M(1,1)*M(2,2)-M(1,2)*M(1,2)
      MI(1,2) = -M(1,2)*M(3,3)+M(1,3)*M(2,3)
      MI(2,1) = MI(1,2)
      MI(1,3) = M(1,2)*M(2,3)-M(1,3)*M(2,2)
      MI(3,1) = MI(1,3)
      MI(2,3) = -M(1,1)*M(2,3)+M(1,3)*M(1,2)
      MI(3,2) = MI(2,3)
      call rscale(9,MI,1./DET)
      return
      END
      
      SUBROUTINE LJGRADPB(Q,GR,N,LJS,LJE,BL)
      IMPLICIT NONE
      INTEGER  I, J, K, N
      real*8   Q(3,N), GR(3,N),  QIJ(3), LJS,  LJE,  RSQ, LJS6,
     &     GRIJ, BL, BL2
      
      LJS6=LJS**6
      call rzero (3*N,GR)
      BL2=BL/2
      
      DO I=1,N-1
         DO J=I+1,N
            RSQ=0
            DO K=1,3
               QIJ(K)=Q(K,I)-Q(K,J)
               IF(QIJ(K) > BL2) QIJ(K)=QIJ(K)-BL
               IF(QIJ(K) < -BL2) QIJ(K)=QIJ(K)+BL
               RSQ=RSQ+QIJ(K)**2
            ENDDO
            DO K=1,3
               GRIJ=4*LJE*LJS6*((12*LJS6*QIJ(K)/RSQ**7)-
     &              (6*QIJ(K)/RSQ**4))
               GR(K,I)=GR(K,I)+GRIJ
               GR(K,J)=GR(K,J)-GRIJ
            ENDDO
         ENDDO
      ENDDO
      return
      END
            


