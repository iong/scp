SUBROUTINE vgw0(Q0, W, TAUMAX,TAUI, Q, G, gamma0)
        use vgw
        use propagation
        use rhss
        IMPLICIT NONE
        REAL*8, intent(in) :: Q0(3,N_atom), TAUMAX, TAUI
        REAL*8, intent(inout) :: Q(3,N_atom), G(3,3,N_atom), gamma0
        REAL*8, intent(out) :: W
        real*8 :: G0(3,3), M(3,3), T, ULJ, LOGDET, DETI, TSTEP=1e-3
        real*8 :: tnow
        integer :: I, nsteps
      

        if (TAUI <= 0.0) then
                T = TAUMIN
                Q=Q0
                call interaction_lists(Q,N_atom,RC,BL, QRC) !Determine which particles

                G0=T*MASS*reshape( (/1,0,0,0,1,0,0,0,1/), (/3, 3/) )

                DO I=1,N_atom
                        G(:,:,I) = G0
                ENDDO

                CALL potential_energy(Q,N_atom,ULJ)
                gamma0=-T*ULJ
        else
                T = 0.5*TAUI
        endif

        call euler(RHSS0, Q, G, gamma0, TSTEP, T, TAUMAX/2.0,1e-4,1e-3)

        LOGDET=0.0
        DO I=1,N_atom
                CALL INVDET(G(:,:,I) , M, DETI)
                LOGDET = LOGDET + LOG(DETI)
        ENDDO
        


        W=-(1/TAUMAX)*(2.0*gamma0 - 0.5*LOGDET)
!      write (*,*) y(1:4)
!      write(*,*) y(2+3*N_atom:7+3*N_atom),  y(2+9*N_atom:10+9*N_atom)
!      write (*,*) y(2+18*N_atom:4+18*N_atom)
        
        
!        the effective potential due to the effective mass should be handled
!        separately
!      W=-(1/TAUMAX)*LNP-((1.5D0/TAUMAX)* &
!              LOG(N_atom*2*3.1415926*TAUMAX/MASS))
        return
END

! vim:et
