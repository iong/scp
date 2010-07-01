SUBROUTINE vgw0(N_ATOM_,IMASS,QCNFG,W, TAUMAX,TAUI,ATOL,RC,Y)
      use vgw
      IMPLICIT NONE
      integer LIW,MF,IERR
      real*8 RTOL,ATOL,ATOMICMASS,RC
      parameter(LIW=20,MF=10,ATOMICMASS=0.020614788876D0, &
           RTOL=0.0D0)
      INTEGER N_atom_,NEQ,LRW,N_STEP,ITASK,ISTATE,IOPT,ITOL,I,J,l, &
           K,CNT,UF,IWORK(LIW)
      REAL*8 QCNFG(3,N_atom_),ENRG,FX(3,N_atom_),QFINAL(3,N_atom_), &
           IMASS,W,LNZ(2),TAU(2),LOGZ,T,C0,UX(3,N_atom_), &
           TOUT,ULJ,LNP,LJS,LJE,BL2, &
           TAUMAX,TSTEP,TAUI,dummy,A(3,3),Z(3,3),F1(3),F2(3),eig(3), &
           Fscaled(3,N_atom_),Fc(3)

      REAL*8 Y(1+9*N_atom_), RWORK(36+336*N_atom_)
      EXTERNAL RHSS0
      EXTERNAL JAC

        N_atom = N_atom_
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
      call interaction_lists(QCNFG,N_atom,RC,BL, QRC) !Determine which particles

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

      CALL potential_energy(QCNFG,N_atom,ULJ)
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


