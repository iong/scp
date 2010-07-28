module propagation
        interface
        SUBROUTINE FF(Q, BLKC, gamma0, Qprime, Gprime, gamma0prime)
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
        interface
        SUBROUTINE F(Q, BLKC, gamma0, Qdot, Gdot, gamma0dot)
                use vgw
                IMPLICIT NONE
                REAL*8, intent(in) :: Q(3,N_atom), BLKC(3,3,N_atom), gamma0
                REAL*8, intent(out) :: Qdot(3,N_atom), Gdot(3,3,N_atom), gamma0dot
        end subroutine
        end interface
 


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

subroutine rk4(F, Q, G, gamma0, dt, nsteps)
        use vgw
        implicit none
        REAL*8, intent(inout) :: Q(3,N_atom), G(3,3,N_atom), gamma0
        real*8, intent(in) :: dt!, tstart, tstop
        integer, intent(in) :: nsteps
        REAL*8 :: Q1(3,N_atom), G1(3,3,N_atom), gamma01
        REAL*8 :: QP(3,N_atom,4), GP(3,3,N_atom,4), gamma0p(4)
        integer :: i
        interface
        SUBROUTINE F(Q, BLKC, gamma0, Qdot, Gdot, gamma0dot)
                use vgw
                IMPLICIT NONE
                REAL*8, intent(in) :: Q(3,N_atom), BLKC(3,3,N_atom), gamma0
                REAL*8, intent(out) :: Qdot(3,N_atom), Gdot(3,3,N_atom), gamma0dot
        end subroutine
        end interface
      

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
end subroutine rk4

subroutine euler(F, Q, G, gamma0, dt, tstart, tstop,atol, rtol)
    use vgw
    implicit none
    REAL*8, intent(inout) :: Q(3,N_atom), G(3,3,N_atom), gamma0
    real*8, intent(in) :: tstart, tstop, atol, rtol
    real*8, intent(inout) :: dt
    REAL*8 :: Q1(3,N_atom), G1(3,3,N_atom), gamma01, rmserr, t
        interface
        SUBROUTINE F(Q, BLKC, gamma0, Qdot, Gdot, gamma0dot)
                use vgw
                IMPLICIT NONE
                REAL*8, intent(in) :: Q(3,N_atom), BLKC(3,3,N_atom), gamma0
                REAL*8, intent(out) :: Qdot(3,N_atom), Gdot(3,3,N_atom), gamma0dot
        end subroutine
        end interface

    t = tstart
    do while (t<tstop)
        call eulerstep(F, Q, G, gamma0, Q1, G1, gamma01, dt, atol, rtol, rmserr)
        !write (*,*) t, dt, rmserr
        if (rmserr <= 1.0) then
            if (t + dt > tstop) then
                exit
            endif
            t = t + dt
            Q=Q1
            G=G1
            gamma0 = gamma01
            if (rmserr <=0.1) then
              dt = dt * 1.52799165273419
            endif
        else
            dt = dt/1.90779652734191
        endif
    enddo
    call eulerstep(F, Q, G, gamma0, Q1, G1, gamma01, tstop-t, atol, rtol, rmserr)
    Q=Q1
    G=G1
    gamma0 = gamma01
    return
end subroutine



subroutine eulerstep(F, Q, G, gamma0, Q1, G1, gamma01, dt, atol, rtol, rmserr)
    use vgw
    implicit none
    REAL*8, intent(in) :: Q(3,N_atom), G(3,3,N_atom), gamma0
    REAL*8, intent(out) :: Q1(3,N_atom), G1(3,3,N_atom), gamma01
    real*8, intent(in) :: dt, atol, rtol
    real*8, intent(out) :: rmserr
    integer, parameter :: p=2
    REAL*8 :: QP(3,N_atom,p), GP(3,3,N_atom,p), gamma0p(p)
    REAL*8 :: Qe(3,N_atom), Ge(3,3,N_atom), gamma0e
    integer :: i, j, k
        interface
        SUBROUTINE F(Q, BLKC, gamma0, Qdot, Gdot, gamma0dot)
                use vgw
                IMPLICIT NONE
                REAL*8, intent(in) :: Q(3,N_atom), BLKC(3,3,N_atom), gamma0
                REAL*8, intent(out) :: Qdot(3,N_atom), Gdot(3,3,N_atom), gamma0dot
        end subroutine
        end interface

    call F(Q, G, gamma0, QP(:,:,1), GP(:,:,:,1), gamma0p(1))
    Q1 = Q + dt*QP(:,:,1)
    G1 = G + dt*GP(:,:,:,1)
    gamma01 = gamma0 + dt*gamma0p(1)
    call F(Q1, G1, gamma01, QP(:,:,2), GP(:,:,:,2), gamma0p(2))

    Q1 = Q + 0.5*dt*(QP(:,:,1) + QP(:,:,2))
    G1 = G + 0.5*dt*(GP(:,:,:,1) + GP(:,:,:,2))
    gamma01 = gamma0 + 0.5*dt*(gamma0P(1) + gamma0P(2))
    Qe = 0.5*dt*(QP(:,:,1) - QP(:,:,2))
    Ge = 0.5*dt*(GP(:,:,:,1) - GP(:,:,:,2))
    gamma0e = 0.5*dt*(gamma0P(1) - gamma0P(2))

    Qe = Qe / (rtol*Q + atol)
    Ge = Ge / (rtol*G + atol)
    gamma0e = gamma0e / (rtol*gamma0 + atol)
    
    rmserr = 0.0
    do i=1,N_atom
        do j=1,3
            do k=1,3
                rmserr = rmserr + Ge(k, j, i)**2
            enddo
            rmserr = rmserr + Qe(j, i)**2
        enddo
    enddo
    rmserr = rmserr! + gamma0e*gamma0e
    rmserr = sqrt(rmserr/(12*N_atom + 1))
    return
   end subroutine

subroutine rk45(F, Q, G, gamma0, dt, tstart, tstop,atol, rtol)
    use vgw
    implicit none
    REAL*8, intent(inout) :: Q(3,N_atom), G(3,3,N_atom), gamma0
    real*8, intent(in) :: tstart, tstop, atol, rtol
    real*8, intent(inout) :: dt
    REAL*8 :: Q1(3,N_atom), G1(3,3,N_atom), gamma01, rmserr, t
        interface
        SUBROUTINE F(Q, BLKC, gamma0, Qdot, Gdot, gamma0dot)
                use vgw
                IMPLICIT NONE
                REAL*8, intent(in) :: Q(3,N_atom), BLKC(3,3,N_atom), gamma0
                REAL*8, intent(out) :: Qdot(3,N_atom), Gdot(3,3,N_atom), gamma0dot
        end subroutine
        end interface

    t = tstart
    do while (t<tstop)
        call rk45step(F, Q, G, gamma0, Q1, G1, gamma01, dt, atol, rtol, rmserr)
        !write (*,*) t, dt, rmserr
        if (rmserr <= 1.0) then
            if (t + dt > tstop) then
                exit
            endif
            t = t + dt
            Q=Q1
            G=G1
            gamma0 = gamma01
            if (rmserr <=0.1) then
                dt = dt * 1.12799165273419
            endif
        else
            dt = dt/1.90779652734191
        endif
    enddo
    call rk45step(F, Q, G, gamma0, Q1, G1, gamma01, tstop-t, atol, rtol, rmserr)
    Q=Q1
    G=G1
    gamma0 = gamma01
    return
end subroutine


subroutine rk45step(F, Q, G, gamma0, Q1, G1, gamma01, dt, atol, rtol, rmserr)
    use vgw
    implicit none
    REAL*8, intent(in) :: Q(3,N_atom), G(3,3,N_atom), gamma0
    REAL*8, intent(out) :: Q1(3,N_atom), G1(3,3,N_atom), gamma01
    real*8, intent(in) :: dt, atol, rtol
    real*8, intent(out) :: rmserr
    integer, parameter :: p=6
    REAL*8 :: QP(3,N_atom,p), GP(3,3,N_atom,p), gamma0p(p)
    REAL*8 :: Qe(3,N_atom), Ge(3,3,N_atom), gamma0e
    integer :: i, j, k
    real*8 :: c(6) = (/ 0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0/)
    real*8 :: b(6) = (/902880.0/7618050.0, 0.0, 3953664.0/7618050.0, &
        3855735.0/7618050.0, -1371249.0/7618050.0, 277020.0/7618050.0 /)
    real*8 :: be(6) = (/ 1.0 / 360.0, 0.0, -128.0 / 4275.0, &
        -2197.0 / 75240.0, 1.0 / 50.0, 2.0 / 55.0 /)
    real*8 :: a(6,6) = reshape( (/ &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        1.0/4.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        3.0/32.0, 9/32.0, 0.0, 0.0, 0.0, 0.0, &
        1932.0/2197.0, -7200.0/2197.0,  7296.0/2197.0, 0.0, 0.0, 0.0, &
        439.0/216.0, -8.0,  3680.0/513.0,    -845.0/4104.0, 0.0, 0.0, &
        -8.0/27.0,  2.0,   -3544.0/2565.0,  1859.0/4104.0, 11.0/40.0, 0.0/), &
            (/ 6, 6 /) )
        interface
        SUBROUTINE F(Q, BLKC, gamma0, Qdot, Gdot, gamma0dot)
                use vgw
                IMPLICIT NONE
                REAL*8, intent(in) :: Q(3,N_atom), BLKC(3,3,N_atom), gamma0
                REAL*8, intent(out) :: Qdot(3,N_atom), Gdot(3,3,N_atom), gamma0dot
        end subroutine
        end interface

    call F(Q, G, gamma0, QP(:,:,1), GP(:,:,:,1), gamma0p(1))
    do j=2,p
        Q1 = Q
        G1=G
        gamma01 = gamma0
        do k=1,j-1
                Q1 = Q1 + a(k,j)*dt*QP(:,:,k)
                G1 = G1 + a(k,j)*dt*GP(:,:,:,k)
                gamma01 = gamma01 + a(j,k)*dt*gamma0p(k)
        enddo
        call F(Q1, G1, gamma01, QP(:,:,j), GP(:,:,:,j), gamma0p(j))
    enddo

    Q1 = Q
    G1=G
    gamma01 = gamma0
    Qe=0.0
    Ge=0.0
    gamma0e=0.0
    do j=1,p
        Q1 = Q1 + dt*b(j)*QP(:,:,j)
        Qe = Qe + dt*be(j)*QP(:,:,j)
        G1 = G1 + dt*b(j)*GP(:,:,:,j)
        Ge = Ge + dt*be(j)*GP(:,:,:,j)
        gamma01 = gamma01 + dt*b(j)*gamma0p(j)
        gamma0e = gamma0e + dt*be(j)*gamma0p(j)
    enddo

    Qe = Qe / (rtol*Q + atol)
    Ge = Ge / (rtol*G + atol)
    gamma0e = gamma0e / (rtol*gamma0 + atol)
    
    rmserr = 0.0
    do i=1,N_atom
        do j=1,3
            do k=1,3
                rmserr = rmserr + Ge(k, j, i)**2
            enddo
            rmserr = rmserr + Qe(j, i)**2
        enddo
    enddo
    rmserr = rmserr! + gamma0e*gamma0e
    rmserr = sqrt(rmserr/(12*N_atom + 1))
    return
   end subroutine rk45step
end module propagation
! vim:et
