SUBROUTINE vgw1(Q0, W, TAUMAX,TAUI, Y, Meff, invMeff)
    use propagation
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(3,N_atom), TAUMAX, TAUI
    REAL*8, intent(inout) :: Y(1+21*N_atom)
    REAL*8, intent(out) :: W, Meff(3,3,N_atom), invMeff(3,3,N_atom)
    real*8 :: G(3,3,N_atom), G0(6), M(3,3), T, ULJ, LOGDET, DETI
    real*8 :: TSTEP=1e-3
    real*8 :: tnow
    integer :: I, nsteps, CNT, CNT2
    external RHSS
  

    if (TAUI <= 0.0) then
        T = TAUMIN

        call interaction_lists(Q0,N_atom,RC,BL, QRC) !Determine which particles
        call Upot_tau0(Q0,N_atom,ULJ)
        call Ux_tau0(Q0, N_atom,Y(2+18*N_atom:1+21*N_atom))
        call init_mylsode(1+21*N_atom)

        Y(1) = -T*ULJ
        Y(2:1+3*N_atom)=reshape(Q0, (/3*N_atom/))
        Y(2+18*N_atom:1+21*N_atom) = Y(2+18*N_atom:1+21*N_atom) * T

        G0=T*invmass*(/1,0,0,1,0,1/)
        CNT=2+3*N_atom
        CNT2=2+9*N_atom
        DO I=1,N_atom
            Y(CNT:CNT+5) = G0
            Y(CNT2:CNT2+8) = (/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/)
            CNT=CNT+6
            CNT2=CNT2+9
        ENDDO
    else
        T = 0.5*TAUI
    endif

    call mylsode(RHSS1, Y, 1+21*N_atom, TSTEP, T, TAUMAX/2.0,atol,rtol)

    call unpack_g(y, G)

    LOGDET=0.0
    DO I=1,N_atom
        CALL INVDET(G(:,:,I) , M, DETI)
        LOGDET = LOGDET + LOG(DETI)
        invMeff(:,:,i) = M/DETI
    ENDDO
    W=-(1/TAUMAX)*(2.0*y(1) - 0.5*LOGDET)

    Meff = G*2.0*mass**2/TAUMAX
    invMeff = invMeff*TAUMAX/(2.0*mass**2)

END SUBROUTINE

! vim:et
