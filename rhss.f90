module rhss
        interface
        SUBROUTINE RHSS0(Q, BLKC, gamma0, Qprime, Gprime, gamma0prime)
              use vgw
              IMPLICIT NONE
              REAL*8, intent(in) :: Q(3,N_atom), BLKC(3,3,N_atom), gamma0
              REAL*8, intent(out) :: Qprime(3,N_atom), Gprime(3,3,N_atom), gamma0prime
        end subroutine
        end interface
end module
