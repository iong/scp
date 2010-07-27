SUBROUTINE vgw0(QCNFG, W, TAUMAX,TAUI,Y)
        use vgw
        IMPLICIT NONE
        integer, parameter :: LIW=20,MF=10
        INTEGER NEQ,LRW,N_STEP,ITASK,ISTATE,IOPT,ITOL,I,J,l, &
                K,CNT,UF,IWORK(LIW)
        REAL*8 QCNFG(3,N_atom),ENRG,FX(3,N_atom),QFINAL(3,N_atom), &
                W,LNZ(2),TAU(2),LOGZ,T,C0,UX(3,N_atom), &
                TOUT,ULJ,LNP,LJS,LJE,BL2, &
                TAUMAX,TSTEP,TAUI
        REAL*8 Y(1+9*N_atom), RWORK(36+336*N_atom)
        integer :: IERR
        real*8, parameter :: RTOL=0.0D0
        EXTERNAL RHSS0
        EXTERNAL JAC

        NEQ=1+9*N_atom
        LRW=36+336*N_atom
        TSTEP=0.1D0
        BL2=BL/2

        
        call DLSODEINIT(IWORK,RWORK,ITASK,IOPT,ISTATE,ITOL,LIW,LRW)
      

        if (T <= 0.0) then
                T = TAUMIN
                call interaction_lists(QCNFG,N_atom,RC,BL, QRC) !Determine which particles
                call requal(3*N_atom,QCNFG,Y(2)) !Copy particle positions to Y

                C0=T*MASS

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

                CALL potential_energy(QCNFG,N_atom,ULJ)
                Y(1)=-T*ULJ
        else
                T = 0.5*TAUI
        endif
        TOUT=TAUMAX/2.0D0


!      write (*,*) y(1:4)
!      write(*,*) y(2+3*N_atom:7+3*N_atom),  y(2+9*N_atom:10+9*N_atom)
!      write (*,*) y(2+18*N_atom:4+18*N_atom)
        
        CALL DLSODE(RHSS0,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE, &
                IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
        CALL LNPS(N_atom,NEQ,Y,LOGZ)

        W=-(1/TAUMAX)*LOGZ
        
!        the effective potential due to the effective mass should be handled
!        separately
!      W=-(1/TAUMAX)*LNP-((1.5D0/TAUMAX)* &
!              LOG(N_atom*2*3.1415926*TAUMAX/MASS))


        return
END

! vim:et
