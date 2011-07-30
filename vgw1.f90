subroutine vgw1(Q0, W, TAUMAX,TAUI, Y, Meff, invMeff)
    use propagation
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(3,N_atom), TAUMAX, TAUI
    REAL*8, intent(inout) :: Y(1+21*N_atom)
    REAL*8, intent(out) :: W
    REAL*8, intent(out), optional:: Meff(3,3,N_atom), invMeff(3,3,N_atom)
    real*8 :: G(3,3,N_atom), M(3,3), T, LOGDET, DETI
    real*8 :: TSTEP=1e-3
    integer :: I
    external RHSS


    if (TAUI <= 0.0d0) then
        call interaction_lists(Q0,N_atom,RC,BL, QRC) !Determine which particles
        call init_gaussians(Q0, TAUMIN, y)
        call init_mylsode(1 + 21*N_atom)

        T = TAUMIN
    else
        T = 0.5*TAUI
    endif

    if (TAUMAX > TAUI) then
        call mylsode(RHSS1, Y, 1+21*N_atom, TSTEP, T, TAUMAX/2.0,atol,rtol)
    end if

    call unpack_g(y, G)

    LOGDET=0.0
    if (present(Meff)) then
        DO I=1,N_atom
            CALL INVDET(G(:,:,I) , M, DETI)
            LOGDET = LOGDET + LOG(DETI)
            invMeff(:,:,i) = M/DETI
        ENDDO
        Meff = G*2.0*mass**2/TAUMAX
        invMeff = invMeff*TAUMAX/(2.0*mass**2)
    else
        DO I=1,N_atom
            CALL INVDET(G(:,:,I) , M, DETI)
            LOGDET = LOGDET + LOG(DETI)
        ENDDO
    end if
    W=-(1/TAUMAX)*(2.0*y(1) - 0.5*LOGDET)

end subroutine vgw1

! vim:et
