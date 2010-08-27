      SUBROUTINE vgwquenchspb(N_ATOM,IMASS,QCNFG,FX,LNP,W,
     &     ENRG,TAUMAX,TAUI,BL,ATOL,RC,Y,
     &     InvMeff,SqrtMeff)
      IMPLICIT NONE
      integer LIW,MF,IERR
      real*8 RTOL,ATOL,ATOMICMASS,RC
      parameter(LIW=20,MF=10,ATOMICMASS=0.020614788876D0,
     &     RTOL=1.0e-4)
      INTEGER N_atom,NEQ,LRW,N_STEP,ITASK,ISTATE,IOPT,ITOL,I,J,l,
     &     K,CNT,UF,IWORK(LIW)
      REAL*8 QCNFG(3,N_atom),ENRG,FX(3,N_atom),QFINAL(3,N_atom),
     &     MASS,IMASS,W,LNZ(2),TAU(2),LOGZ,T,C0,UX(3,N_atom),
     &     TOUT,ULJ,LNP,LJS,LJE,BL,BL2,
     &     TAUMAX,TSTEP,TAUI,dummy,A(3,3),Z(3,3),F1(3),F2(3),eig(3),
     &     Fscaled(3,N_atom),Fc(3)

      REAL*8 Y(1+21*N_atom), RWORK(36+336*N_atom),SqrtMeff(3,3,N_ATOM),
     &     InvMeff(3,3,N_ATOM),meff,SqrtMeffInv(3,3,N_atom)
      EXTERNAL RHSS
      EXTERNAL JAC
      
      NEQ=1+21*N_atom
      LRW=36+336*N_atom
      TSTEP=0.1D0
      BL2=BL/2
                
      MASS=1.0D0/(ATOMICMASS*IMASS)

      call requal(3*N_atom,QCNFG,Y(2)) !Copy particle positions to Y
      
      call initialize_SG_3G(MASS,BL,N_atom) 
      call initialize_QRC(QCNFG,N_atom,RC,BL) !Determine which particles interact having folded particles into a box of side bl
 
      CALL DLSODEINIT(IWORK,RWORK,ITASK,IOPT,ISTATE,ITOL,LIW,LRW)

      C0=TAUI*MASS
      
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
      W=-(1/TAUMAX)*LNP
    
      call requal(3*N_atom,Y(2+18*N_atom),FX)
      call rscale(3*N_atom,FX,1/TAU(2))
              
! effective mass
      dummy=2*IMASS**2/taumax*ATOMICMASS
      meff=0d0
      CNT=2+3*N_atom
      call rzero(3,Fc)
      call rzero(9*N_atom,SqrtMeffInv)
      call rzero(9*N_atom,SqrtMeff)
      call rzero(9*N_atom,InvMeff)
      DO l=1,N_atom            
         DO i=1,3
            DO j=i,3
               A(i,j)=Y(CNT)*dummy
               A(j,i)=A(i,j)
               CNT=CNT+1
            ENDDO
         ENDDO
         call RS(3,3,A,eig,1,Z,F1,F2,IERR)
c         write(6,*) A
c         write(6,*)
c         write(6,*) eig
c         write(6,*)
c         write(6,*) z
c         write(6,*)
         if(IERR.ne.0) stop 'IERR ne 0'
         do j=1,3
            meff=meff+eig(j)
            do i=1,3
               do k=1,3
                  InvMeff(i,j,l)=InvMeff(i,j,l)
     &                 +Z(i,k)*Z(j,k)/eig(k)
               enddo
            enddo
         enddo
c         write(6,*) 'l=',l
c         write(6,*) ((InvMeff(i,j,l), i=1,3),j=1,3)
         do k=1,3
            eig(k)=dsqrt(eig(k))/sqrt(IMASS)
         enddo
         do j=1,3
            do i=1,3
               do k=1,3
                  SqrtMeff(i,j,l)=SqrtMeff(i,j,l)
     &                 +eig(k)*Z(i,k)*Z(j,k)
                  SqrtMeffInv(i,j,l)=SqrtMeffInv(i,j,l)
     &                 +Z(i,k)*Z(j,k)/eig(k)
               enddo
            enddo
         enddo
      enddo

      meff=meff/(N_atom*3)
      
 
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
      INTEGER NGAUSS,N_atom,I,J,K,I1,I2,IG,NEQ,CNT,CNT2,MAXN
      PARAMETER (NGAUSS=4,MAXN=5000)
      REAL*8 T,AG(3,3),GU(3,3),Y(NEQ+3*(NEQ-1)/21),YPRIME(NEQ),
     &     DETA,DETAG,GUG,QP,TRUXXGI,FACTOR,U,UX,UXX,QZQ,EXPAV,
     &     TRMG,DETS,DETI,MASS,LJS,LJE,GUQ,
     &     BL,BL2,M(3,3),A(3,3), 
     &     R(3), Z(3,3),Q12(3)
      REAL*8 Q(3,(NEQ-1)/21),UPV(3,(NEQ-1)/21),BLKC(3,3,(NEQ-1)/21),
     &     UPM(3,3,(NEQ-1)/21), QNK(3,3,(NEQ-1)/21) 
      
      real*8 LJC(NGAUSS),LJA(NGAUSS)
      logical QRC(MAXN*MAXN)
      common /LJ3/ LJC,LJA
      common /systemee/ MASS,BL,N_atom
      common /carray/ QRC
c     equivalence (Y(1),gamma), (Y(2),Q(1,1))

      if(N_atom.ne.(NEQ-1)/21) then
        write (*,*) N_atom, (NEQ-1)/21
        stop 'RHSS: Oops!'
      endif
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

      subroutine initialize_QRC(Q,N,rc,bl)
      implicit none
      integer N,I,J,K,MAXN
      real*8 Q(3,N),rsq,rc,rc2,bl,bl2,qij
      parameter(MAXN=5000)
      logical QRC(MAXN*MAXN)
      common /carray/ QRC
      if(MAXN<N) stop 'initialize_QRC: increase MAXN'

      rc2=rc**2
      bl2=bl/2

      do I=1,N-1
        do J=I+1,N
          rsq=0.0D0
          do K=1,3
            qij=Q(K,I)-Q(K,J)
            if(qij > bl2) qij=qij-bl
            if(qij < -bl2) qij=qij+bl
            rsq=rsq+qij**2
          enddo
          if(rsq.GE.rc2) then
            QRC(I+(J-1)*N)=.FALSE.
          else
            QRC(I+(J-1)*N)=.TRUE.
          endif
          QRC(J+(I-1)*N)=QRC(I+(J-1)*N)
        enddo
      enddo
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
              GRAD=GRAD+2.0*LJC(P)*LJA(P)*EXP(-RSQ*LJA(P))
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
         A(i)=0.0
      enddo
      return
      end

      subroutine matvec(N,A,B,C)
      implicit none
      integer i,N,j
      real*8 A(N,N),B(N),C(N)
      do i=1,N
         C(i)=0
         do j=1,N
            C(i)=C(i)+A(j,i)*B(j)
         enddo
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
            
      
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
C
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,NM,IERR,MATZ
      DOUBLE PRECISION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        A  CONTAINS THE REAL SYMMETRIC MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)
      DOUBLE PRECISION F,G,H,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
  100 CONTINUE
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(D(K))
C
         IF (SCALE .NE. 0.0D0) GO TO 140
C
         DO 125 J = 1, L
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = 0.0D0
  125    CONTINUE
C
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 300
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         IF (L .EQ. 1) GO TO 285
C     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0D0
C
         DO 240 J = 1, L
            F = D(J)
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + A(K,J) * D(K)
               E(K) = E(K) + A(K,J) * F
  200       CONTINUE
C
  220       E(J) = G
  240    CONTINUE
C     .......... FORM P ..........
         F = 0.0D0
C
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
C
         H = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - H * D(J)
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
C
            DO 260 K = J, L
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)
C
  280    CONTINUE
C
  285    DO 290 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F * SCALE
  290    CONTINUE
C
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
      DOUBLE PRECISION F,G,H,HH,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
C          PRODUCED IN THE REDUCTION.
C
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
C
         DO 80 J = I, N
   80    Z(J,I) = A(J,I)
C
         D(I) = A(N,I)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 510
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(D(K))
C
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = D(L)
C
         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  135    CONTINUE
C
         GO TO 290
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         F = D(L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
C     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0D0
C
         DO 240 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + Z(K,J) * D(K)
               E(K) = E(K) + Z(K,J) * F
  200       CONTINUE
C
  220       E(J) = G
  240    CONTINUE
C     .......... FORM P ..........
         F = 0.0D0
C
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
C
            DO 260 K = J, L
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
C
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  280    CONTINUE
C
  290    D(I) = H
  300 CONTINUE
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         H = D(I)
         IF (H .EQ. 0.0D0) GO TO 380
C
         DO 330 K = 1, L
  330    D(K) = Z(K,I) / H
C
         DO 360 J = 1, L
            G = 0.0D0
C
            DO 340 K = 1, L
  340       G = G + Z(K,I) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  360    CONTINUE
C
  380    DO 400 K = 1, L
  400    Z(K,I) = 0.0D0
C
  500 CONTINUE
C
  510 DO 520 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0D0
  520 CONTINUE
C
      Z(N,N) = 1.0D0
      E(1) = 0.0D0
      RETURN
      END
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      DOUBLE PRECISION D(N),E2(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E2 HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = PYTHAG(P,1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0D0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TRBAK1(NM,N,A,E,M,Z)
C
      INTEGER I,J,K,L,M,N,NM
      DOUBLE PRECISION A(NM,N),E(N),Z(NM,M)
      DOUBLE PRECISION S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED1.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  TRED1
C          IN ITS STRICT LOWER TRIANGLE.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK1 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
         L = I - 1
         IF (E(I) .EQ. 0.0D0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0D0
C
            DO 110 K = 1, L
  110       S = S + A(I,K) * Z(K,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / A(I,L)) / E(I)
C
            DO 120 K = 1, L
  120       Z(K,J) = Z(K,J) + S * A(I,K)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END

      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D0
      TST1 = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
         H = DABS(D(L)) + DABS(E(L))
         IF (TST1 .LT. H) TST1 = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + DABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + DSIGN(R,P))
         D(L1) = E(L) * (P + DSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + DABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END

      DOUBLE PRECISION FUNCTION EPSLON (X)
      DOUBLE PRECISION X
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
      DOUBLE PRECISION A,B,C,EPS
C
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT
C            NUMBERS IS NOT A POWER OF THREE.
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C            THE ACCURACY USED IN FLOATING POINT VARIABLES
C            THAT ARE STORED IN MEMORY.
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C     ASSUMPTION 2.
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C            C  IS NOT EXACTLY EQUAL TO ONE,
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM
C                 THE NEXT LARGER FLOATING POINT NUMBER.
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
C
C     THIS VERSION DATED 4/6/83.
C
      A = 4.0D0/3.0D0
   10 B = A - 1.0D0
      C = B + B + B
      EPS = DABS(C-1.0D0)
      IF (EPS .EQ. 0.0D0) GO TO 10
      EPSLON = EPS*DABS(X)
      RETURN
      END

