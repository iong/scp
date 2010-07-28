module propagation
contains
subroutine verlet(F, Q, G, gamma0, tstart, tstop, dt)
        use vgw
      REAL*8, intent(inout) :: Q(3,N_atom), G(3,3,N_atom), gamma0
        real*8, intent(in) :: tstart, tstop, dt
      REAL*8 :: Q1(3,N_atom), G1(3,3,N_atom)
      REAL*8 :: QP(3,N_atom), GP(3,3,N_atom), gamma0p

        nsteps = (tstop - tstart) / dt
        do i = 1, nsteps
                call F(Q, G, gamma0, QP, GP, gamma0p)
                Q1 = Q + 0.5*dt*QP
                G1 = G + 0.5*dt*GP
                call F(N_atom, Q1, G1, gamma0, QP, GP, gamma0p)
                Q = Q + dt*QP
                G = G + dt*GP
                gamma0 = gamma0 + dt*gamma0p
        enddo
end subroutine verlet
end module propagation
! vim:et
