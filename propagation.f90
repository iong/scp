module propagation
        interface
        SUBROUTINE F(Q, BLKC, gamma0, Qprime, Gprime, gamma0prime)
                use vgw
              IMPLICIT NONE
              REAL*8, intent(in) :: Q(3,N_atom), BLKC(3,3,N_atom), gamma0
              REAL*8, intent(out) :: Qprime(3,N_atom), Gprime(3,3,N_atom), gamma0prime
        end subroutine
        end interface
 
contains
subroutine verlet(F, Q, G, gamma0, dt, nsteps)
        use vgw
        implicit none
        REAL*8, intent(inout) :: Q(3,N_atom), G(3,3,N_atom), gamma0
        real*8, intent(in) :: dt!, tstart, tstop
        integer, intent(in) :: nsteps
        REAL*8 :: Q1(3,N_atom), G1(3,3,N_atom)
        REAL*8 :: QP(3,N_atom), GP(3,3,N_atom), gamma0p
        integer :: i

        do i = 1, nsteps
                call F(Q, G, gamma0, QP, GP, gamma0p)
                Q1 = Q + 0.5*dt*QP
                G1 = G + 0.5*dt*GP
                call F(Q1, G1, gamma0, QP, GP, gamma0p)
                Q = Q + dt*QP
                G = G + dt*GP
                gamma0 = gamma0 + dt*gamma0p
        enddo
end subroutine verlet

subroutine rk45(F, Q, G, gamma0, dt, nsteps)
        use vgw
        implicit none
        REAL*8, intent(inout) :: Q(3,N_atom), G(3,3,N_atom), gamma0
        real*8, intent(in) :: dt!, tstart, tstop
        integer, intent(in) :: nsteps
        REAL*8 :: Q1(3,N_atom), G1(3,3,N_atom), gamma01
        REAL*8 :: QP(3,N_atom,4), GP(3,3,N_atom,4), gamma0p(4)
        integer :: i
      

        do i = 1, nsteps
                call F(Q, G, gamma0, QP(:,:,1), GP(:,:,:,1), gamma0p(1))

                Q1 = Q + 0.5*dt*QP(:,:,1)
                G1 = G + 0.5*dt*GP(:,:,:,1)
                gamma01 = gamma0 + 0.5*dt*gamma0p(1)
                call F(Q1, G1, gamma01, QP(:,:,2), GP(:,:,:,2), gamma0p(2))

                Q1 = Q + 0.5*dt*QP(:,:,2)
                G1 = G + 0.5*dt*GP(:,:,:,2)
                gamma01 = gamma0 + 0.5*dt*gamma0p(2)
                call F(Q1, G1, gamma01, QP(:,:,3), GP(:,:,:,3), gamma0p(3))

                Q1 = Q + dt*QP(:,:,3)
                G1 = G + dt*GP(:,:,:,3)
                gamma01 = gamma0 + dt*gamma0p(3)
                call F(Q1, G1, gamma01, QP(:,:,4), GP(:,:,:,4), gamma0p(4))

                Q = Q + (dt/6.0)*(QP(:,:,1) + 2.0*(QP(:,:,2) + QP(:,:,3)) + QP(:,:,4))
                G = G + (dt/6.0)*(GP(:,:,:,1) + 2.0*(GP(:,:,:,2) + GP(:,:,:,3)) + GP(:,:,:,4))
                gamma0 = gamma0 + (dt/6.0)*(gamma0p(1) + 2.0*(gamma0p(2) + gamma0p(3)) + gamma0p(4))
        enddo
end subroutine rk45


subroutine rk(F, Q, G, gamma0, dt, nsteps)
        use vgw
        implicit none
        REAL*8, intent(inout) :: Q(3,N_atom), G(3,3,N_atom), gamma0
        real*8, intent(in) :: dt!, tstart, tstop
        integer, intent(in) :: nsteps
        REAL*8 :: Q1(3,N_atom), G1(3,3,N_atom), gamma01
        REAL*8 :: QP(3,N_atom,4), GP(3,3,N_atom,4), gamma0p(4)
        integer :: i
      

        do i = 1, nsteps
            call F(Q, G, gamma0, QP(:,:,1), GP(:,:,:,1), gamma0p(1))
            do j=2,p
                Q1 = Q
                G1=G
                gamma01 = gamma0
                do k=1,j-1
                        Q1 = Q1 + a(j,k)*dt*QP(:,:,k)
                        G1 = G1 + a(j,k)*dt*GP(:,:,:,k)
                        gamma01 = gamma01 + a(j,k)*dt*gamma0p(k)
                enddo
                call F(Q1, G1, gamma01, QP(:,:,j), GP(:,:,:,j), gamma0p(j))
            enddo

            do j=1,p
                Q = Q + dt*b(j)*QP(:,:,j)
                G = G + dt*b(j)*GP(:,:,:,j)
                G = G + dt*b(j)*gamma0p(j)
            enddo
        enddo
end subroutine rk45




end module propagation
! vim:et
