SUBROUTINE vgw0(Q0, W, TAUMAX,TAUI, Y)
    use propagation
    use utils
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(3,N_atom), TAUMAX(:), TAUI
    REAL*8, intent(inout) :: Y(1+9*N_atom)
    REAL*8, intent(out) :: W(:)
    real*8 :: G(3,3,N_atom), G0(6), M(3,3), T, ULJ, LOGDET, DETI
    real*8 :: TSTEP=1e-3
    real*8 :: tnow
    integer :: I, j, nsteps, CNT, CNT2
  

    if (TAUI <= 0.0) then
        T = TAUMIN

        call interaction_lists(Q0,N_atom,RC,BL, QRC) !Determine which particles
        call Upot_tau0(Q0,N_atom,ULJ)
        call init_mylsode(1+9*N_atom)

        Y(1) = -T*ULJ
        Y(2:1+3*N_atom)=reshape(Q0, (/3*N_atom/))

        G0=T*invmass*(/1,0,0,1,0,1/)
        CNT=2+3*N_atom
        DO I=1,N_atom
            Y(CNT:CNT+5) = G0
            CNT=CNT+6
        ENDDO

    else
        T = 0.5*TAUI
    endif

    do i=1,size(TAUMAX)
        call mylsode(RHSS0, Y, 1+9*N_atom, TSTEP, T, TAUMAX(i)/2.0,atol,rtol)

        call unpackg(N_atom, y(2+3*N_atom:1+9*N_atom), G)
        LOGDET=0.0
        DO j=1,N_atom
            CALL INVDET(G(:,:,j) , M, DETI)
            LOGDET = LOGDET + LOG(DETI)
        ENDDO
    


        W(i)=-(1/TAUMAX(i))*(2.0*y(1) - 0.5*LOGDET - 3.0*N_atom*log(2.0*sqrt(M_PI)))
    enddo
END SUBROUTINE


SUBROUTINE vgw0s(Q0, W, TAUMAX,TAUI, Y)
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(3,N_atom), TAUMAX, TAUI
    REAL*8, intent(inout) :: Y(1+9*N_atom)
    REAL*8, intent(out) :: W
    real*8, dimension(1) :: W_, TAUMAX_

    TAUMAX_(1) = TAUMAX
    call vgw0(Q0, W_, TAUMAX_, TAUI, Y)
    W = W_(1)
END SUBROUTINE

! vim:et
