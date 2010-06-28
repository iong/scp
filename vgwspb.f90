module vgw
        implicit none
        save
        integer :: NGAUSS, N_atom, NEQ
        real*8 , allocatable :: LJC(:),LJA(:)
        real*8 :: MASS, BL
        logical, allocatable :: QRC(:)
contains

SUBROUTINE JAC 
      return
END SUBROUTINE
      

real*8 function DET3(A) result(det)
        real*8 :: A(3,3)
        det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
                + A(2,1)*A(3,2)*A(1,3) - A(1,3)*A(2,2)*A(3,1) &
                - A(1,2)*A(2,1)*A(3,3) - A(2,3)*A(3,2)*A(1,1)
end function


SUBROUTINE lnrho(N_atom,NEQ,Y,LOGZ)
      implicit none
      integer :: I,J,K,N_atom,NEQ,CNT
      real*8 :: Y(NEQ), C(3,3), LOGZ, LOGDET
      
      CNT=2+3*N_atom
      LOGDET=0.0D0
      
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
        LOGDET = LOGDET + LOG(DET3(C))
      ENDDO
      
      LOGZ = 2.0D0*Y(1) - 0.5*LOGDET
      return
END SUBROUTINE
        
      
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
END SUBROUTINE
      
SUBROUTINE RHSS(NEQ,T,Y,YPRIME)
        IMPLICIT NONE
        integer :: I,J,K,I1,I2,IG,NEQ,CNT,CNT2
        REAL*8 :: T,AG(3,3),GU(3,3),Y(NEQ+3*N_atom),YPRIME(NEQ), INVMASS, &
                DETA, DETAG, &
                GUG,QP,TRUXXGI,FACTOR,U,UX,UXX,QZQ,EXPAV, TRMG,DETS,DETI,GUQ, &
                BL2,M(3,3),A(3,3), R(3), Z(3,3),Q12(3)
        REAL*8 Q(3,N_atom),UPV(3,N_atom),BLKC(3,3,N_atom), UPM(3,3,N_atom), &
                QNK(3,3,N_atom) 
     

        BL2=BL/2
        INVMASS = 1.0/MASS
     
        Q = reshape(Y(2:1+3*N_atom), (/3, N_atom/))
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
        A = BLKC(1:3,1:3,I)
        CALL INVDET(A,M,DETS)
        DETS=1.0D0/DETS
         
        DO J=1,3
          TRMG=TRMG + M(J,J)*DETS
        ENDDO
      ENDDO
      
      TRMG = TRMG*INVMASS
      
        U=0.0D0
        UPV = 0.0
        UPM = 0.0
      
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
      
       YPRIME(CNT) = -0.250*TRUXXGI-U
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
              GUG=GUG+INVMASS
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
      
      SUBROUTINE CLRADIUS(Q, N, RAD)
      IMPLICIT NONE
      INTEGER :: N, I, J
      REAL*8 :: Q(3,N), RCM(3), CV(3), CR, RAD
      
      RAD=0
        RCM = 0.0
      
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
END SUBROUTINE


subroutine initialize_QRC(Q, N, RC)
      implicit none
      integer :: N,I,J,K
      real*8 :: Q(3,N),rsq,rc,rc2,bl,bl2,qij

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
end SUBROUTINE

SUBROUTINE GAUSSENERGYPB(Q,N,BL,U,UX)
        IMPLICIT NONE
        INTEGER, intent(in) :: N
        REAL*8, intent(in) :: Q(3,N), BL
        real*8 ,intent(out) :: U,Ux(3, N)
        INTEGER :: I,J,K,P
        real*8 :: RSQ,BL2,QIJ,VEC(3),GRAD
      
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

        Ux = 0.0
          
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
END SUBROUTINE
      
SUBROUTINE FORCES(N_atom,NEQ,Y,FX,Q,TAU)
      IMPLICIT NONE
      INTEGER :: I,J,K,IND,N_atom,NEQ
      REAL*8  :: TAU,Y(NEQ),FX(3,N_atom),Q(3,N_atom),QX(3,N_atom), BLK(3,3),INV(3,3)
      
      
        BLK=0.0
        FX=0.0
      
        IND=2
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
END SUBROUTINE

SUBROUTINE INVS(M,MI)
      IMPLICIT NONE
      REAL*8 :: DET,M(3,3),MI(3,3)
      
      DET=M(1,1)*(M(2,2)*M(3,3)-M(3,2)*M(2,3))-M(1,2)*(M(2,1)*M(3,3)- &
          M(3,1)*M(2,3))+M(1,3)*(M(1,2)*M(3,2)-M(3,1)*M(2,2))
      
      MI(1,1) = M(2,2)*M(3,3)-M(2,3)*M(2,3)
      MI(2,2) = M(1,1)*M(3,3)-M(1,3)*M(1,3)
      MI(3,3) = M(1,1)*M(2,2)-M(1,2)*M(1,2)
      MI(1,2) = -M(1,2)*M(3,3)+M(1,3)*M(2,3)
      MI(2,1) = MI(1,2)
      MI(1,3) = M(1,2)*M(2,3)-M(1,3)*M(2,2)
      MI(3,1) = MI(1,3)
      MI(2,3) = -M(1,1)*M(2,3)+M(1,3)*M(1,2)
      MI(3,2) = MI(2,3)

        MI = MI/DET
      return
END SUBROUTINE
      
SUBROUTINE LJGRADPB(Q,GR,N,LJS,LJE,BL)
      IMPLICIT NONE
      INTEGER  :: I, J, K, N
      real*8 ::  Q(3,N), GR(3,N),  QIJ(3), LJS,  LJE,  RSQ, LJS6, &
          GRIJ, BL, BL2
      
      LJS6=LJS**6
        GR = 0.0
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
               GRIJ=4*LJE*LJS6*((12*LJS6*QIJ(K)/RSQ**7)-(6*QIJ(K)/RSQ**4))
               GR(K,I)=GR(K,I)+GRIJ
               GR(K,J)=GR(K,J)-GRIJ
            ENDDO
         ENDDO
      ENDDO
      return
END SUBROUTINE
end module vgw

! LJC = (/96609.488289873d0, 14584.62075507514d0, -365.460614956589d0, -19.5534697800036d0/)
! LJA = (/1.038252215127D0, 0.5974039109464D0, 0.196476572277834D0, 0.06668611771781D0/)
subroutine vgwinit(N_atom_, ngauss_, mass_, lja_, ljc_, bl_)
        use vgw
        implicit none
        integer, intent(in) :: N_atom_, ngauss_
        real*8, intent(in) :: mass_, bl_
        real*8, intent(in), dimension(ngauss_) :: lja_, ljc_

        N_atom = N_atom_
        NGAUSS = ngauss_
        MASS = mass_
        BL = bl_
        allocate (LJA(NGAUSS), LJC(NGAUSS), QRC(N_atom*N_atom))
        LJA(1:NGAUSS) = lja_(1:ngauss)
        LJC(1:NGAUSS) = ljc_(1:ngauss)
        
        neq = 1+21*N_atom
end subroutine


subroutine vgwquenchspb(QCNFG,FX, W, ENRG,TAUMAX,TAUI,ATOL,RC,Y)
        use vgw
        IMPLICIT NONE
        integer, parameter :: LIW=20, MF=10
        real*8 ::  ATOL,RC
        INTEGER :: LRW,ITASK,ISTATE,IOPT,ITOL,I,CNT
        REAL*8 :: QCNFG(3,N_atom),ENRG,FX(3,N_atom), W,LNZ(2),TAU(2),LOGZ,T,&
                        TOUT,ULJ,BL2, TAUMAX,TSTEP,TAUI
        REAL*8 :: Y(NEQ), RWORK(36+336*N_atom)
        integer :: IWORK(LIW)
        integer :: ANEQ(2)
        real*8 :: RTOL(2)

    

      LRW=36+336*N_atom
      TSTEP=0.1D0
      BL2=BL/2
      
             !Determine which particles interact having folded particles into a box of side bl
        call initialize_QRC(QCNFG, N_atom, RC)
 
        ANEQ = NEQ
        ITOL=1
        RTOL=0.0
        ITASK=1
        IOPT=1
        ISTATE=1
        IWORK(5:10)=(/4, 100000, 0, 0, 0, 0/)
        RWORK(5:10)=0.0D0

        ! INITIALIZE Y VECTOR WITH INITITIAL CONDITIONS
        Y(2:1+3*N_atom) = reshape(QCNFG, (/3*N_atom/))
        CNT=2+3*N_atom
        DO I=1,N_atom             
                Y(CNT:CNT+5)=(/1.0, 0.0, 0.0, 1.0, 0.0, 1.0/)*TAUI/MASS
                cnt = cnt + 6
        ENDDO
      
        DO I=1,N_atom             ! INITIALIZE Q MATRIX
                Y(CNT:CNT+8)=(/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/)
                cnt = cnt + 9
        ENDDO
    
        CALL GAUSSENERGYPB(QCNFG,N_atom,BL,ULJ,Y(CNT))      
        y(cnt:cnt+3*N_atom-1)=y(cnt:cnt+3*N_atom-1)*taui
        Y(1)=-TAUI*ULJ

        T=TAUI
        IF(TAUMAX < 1) TSTEP = 0.001D0
        IF((TAUMAX/2)-TAUI < TSTEP) TSTEP = ((TAUMAX/2)-TAUI)/2.0D0
     
        TAU(1) = TAUMAX/2.0D0 - TSTEP
        TAU(2) = TAUMAX/2.0D0

     
        DO I=1,2
                TOUT=TAU(I)
                CALL DLSODE(RHSS,ANEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE, &
                        IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
                CALL lnrho(N_atom,NEQ,Y,LOGZ)
                LNZ(I) = LOGZ+1.5D0*LOG(TOUT)
        ENDDO
  
        ENRG=-(LNZ(2)-LNZ(1))/(2*TSTEP)
        W=-(LOGZ - 1.5D0 * N_atom * LOG(2*3.1415926*TAUMAX/MASS) ) / TAUMAX
    
        fx = reshape(y(2+18*N_atom : 1 + 21*N_atom), (/3, N_atom/)) / TAU(2)

 
      return
END SUBROUTINE

subroutine effective_mass(beta, y, InvMeff, SqrtMeff)
        use vgw
        implicit none
        real*8, intent(in) :: beta, y(neq)
        real*8, intent(out) :: InvMeff(3,3,N_atom), SqrtMeff(3,3,N_ATOM)
        real*8 :: meff, SqrtMeffInv(3,3,N_atom), msqkbt,A(3,3),Z(3,3),F1(3),F2(3),eig(3)
        integer :: CNT, i, j, k, l, IERR
              
! effective mass
      msqkbt=2*MASS**2/beta
      meff=0d0
      CNT=2+3*N_atom

        SqrtMeffInv = 0.0
        SqrtMeff = 0.0
        InvMeff = 0.0
      DO l=1,N_atom            
         DO i=1,3
            DO j=i,3
               A(i,j)=Y(CNT)*msqkbt
               A(j,i)=A(i,j)
               CNT=CNT+1
            ENDDO
         ENDDO
                call RS(3,3,A,eig,1,Z,F1,F2,IERR)
                if(IERR.ne.0) stop 'IERR ne 0'
         do j=1,3
            meff=meff+eig(j)
            do i=1,3
               do k=1,3
                  InvMeff(i,j,l)=InvMeff(i,j,l) + Z(i,k)*Z(j,k)/eig(k)
               enddo
            enddo
         enddo
         do k=1,3
            eig(k)=dsqrt(eig(k))
         enddo
         do j=1,3
            do i=1,3
               do k=1,3
                  SqrtMeff(i,j,l)=SqrtMeff(i,j,l) + eig(k)*Z(i,k)*Z(j,k)
                  SqrtMeffInv(i,j,l)=SqrtMeffInv(i,j,l) + Z(i,k)*Z(j,k)/eig(k)
               enddo
            enddo
         enddo
      enddo

      meff=meff/(N_atom*3)
end subroutine
 
