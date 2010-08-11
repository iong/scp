module propagation
contains
subroutine verlet(F, x, NEQ, dt, nsteps)
    implicit none
    integer, intent(in) :: NEQ
    REAL*8, intent(inout) :: x(NEQ)
    real*8, intent(in) :: dt!, tstart, tstop
    integer, intent(in) :: nsteps
    REAL*8 :: x1(NEQ), xp(NEQ)
    integer :: i
    external F

    do i = 1, nsteps
            call F(x, xp)
            x1 = x + 0.5*dt*xp
            call F(x1, xp)
            x = x + dt*xp
    enddo
end subroutine verlet

subroutine rk4(F, x, NEQ, dt, nsteps)
    implicit none
    integer, intent(in) :: NEQ
    REAL*8, intent(inout) :: x(NEQ)
    real*8, intent(in) :: dt!, tstart, tstop
    integer, intent(in) :: nsteps
    REAL*8 :: x1(NEQ), xp(NEQ,5)
    integer :: i
  

    do i = 1, nsteps
            call F(x, xp(:,1))

            x1 = x + 0.5*dt*xp(:,1)
            call F(x1, xp(:,2))

            x1 = x + 0.5*dt*xp(:,2)
            call F(x1, xp(:,3))

            x1 = x + dt*xp(:,3)
            call F(x1, xp(:,4))

            x = x + (dt/6.0)*(xp(:,1) + 2.0*(xp(:,2) + xp(:,3)) + xp(:,4))
    enddo
end subroutine rk4

subroutine euler(F, x, NEQ, dt, tstart, tstop,atol, rtol)
    implicit none
    integer, intent(in) :: NEQ
    REAL*8, intent(inout) :: x(NEQ)
    real*8, intent(in) :: tstart, tstop, atol, rtol
    real*8, intent(inout) :: dt
    REAL*8 :: x1(NEQ), rmserr, t
    external F

    t = tstart
    do while (t<tstop)
        call eulerstep(F, x, x1, NEQ, dt, atol, rtol, rmserr)
        !write (*,*) t, dt, rmserr
        if (rmserr <= 1.0) then
            if (t + dt > tstop) then
                exit
            endif
            t = t + dt
            x=x1
            if (rmserr <=0.1) then
              dt = dt * 1.52799165273419
            endif
        else
            dt = dt/1.90779652734191
        endif
    enddo
    call eulerstep(F, x, x1, NEQ, tstop-t, atol, rtol, rmserr)
    x=x1
    return
end subroutine



subroutine eulerstep(F, x, x1, NEQ, dt, atol, rtol, rmserr)
    implicit none
    integer, intent(in) :: NEQ
    REAL*8, intent(in) :: x(NEQ)
    REAL*8, intent(out) :: x1(NEQ)
    real*8, intent(in) :: dt, atol, rtol
    real*8, intent(out) :: rmserr
    integer, parameter :: p=2
    REAL*8 :: xp(NEQ,p), xe(NEQ)
    integer :: i, j, k
    external F

    call F(x, xp(:,1))
    x1 = x + dt*xp(:,1)
    call F(x1, xp(:,2))

    x1 = x + 0.5*dt*(xp(:,1) + xp(:,2))
    xe = 0.5*dt*(xp(:,1) - xp(:,2))

    xe = xe / (abs(rtol*x) + atol)
    rmserr = sqrt( sum(xe**2) / NEQ )
    
    return
   end subroutine

subroutine rk45(F, x, NEQ, dt, tstart, tstop,atol, rtol)
    implicit none
    integer, intent(in) :: NEQ
    REAL*8, intent(inout) :: x(NEQ)
    real*8, intent(in) :: tstart, tstop, atol, rtol
    real*8, intent(inout) :: dt
    REAL*8 :: x1(NEQ), rmserr, t
    external F

    t = tstart
    do while (t<tstop)
        call rk45step(F, x, x1, NEQ, dt, atol, rtol, rmserr)
        !write (*,*) t, dt, rmserr
        if (rmserr <= 1.0) then
            if (t + dt > tstop) then
                exit
            endif
            t = t + dt
            x=x1
            if (rmserr <=0.1) then
              dt = dt * 1.52799165273419
            endif
        else
            dt = dt/1.90779652734191
        endif
    enddo
    call rk45step(F, x, x1, NEQ, tstop-t, atol, rtol, rmserr)
    x=x1
    return
end subroutine


subroutine rk45step(F, x, x1, NEQ, dt, atol, rtol, rmserr)
    implicit none
    integer, intent(in) :: NEQ
    REAL*8, intent(in) :: x(NEQ)
    REAL*8, intent(out) :: x1(NEQ)
    real*8, intent(in) :: dt, atol, rtol
    real*8, intent(out) :: rmserr
    integer, parameter :: p=6
    REAL*8 :: xp(NEQ,p), xe(NEQ)
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
    external F

    call F(x, xp(:,1))
    do j=2,p
        x1 = x
        do k=1,j-1
                x1 = x1 + a(k,j)*dt*xp(:,k)
        enddo
        call F(x1, xp(:,j))
    enddo

    x1 = x
    xe=0.0
    do j=1,p
        x1 = x1 + dt*b(j)*xp(:,j)
        xe = xe + dt*be(j)*xp(:,j)
    enddo
    xe = xe / (abs(rtol*x) + atol)
    rmserr = sqrt( sum(xe**2) / NEQ )

    return
   end subroutine rk45step
end module propagation
! vim:et:softtabstop=4:sw=4
