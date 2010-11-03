module propagation
    integer, parameter :: LIW=20, MF=10
    integer :: ITASK,IOPT,ISTATE,LRW,ITOL,IWORK(LIW),IERR
    real*8, allocatable :: RWORK(:)


contains
subroutine  init_mylsode(NEQ)
    integer, intent(IN) :: NEQ
    LRW = 10 + 16*NEQ
    if (.not. allocated(RWORK)) then
        allocate(RWORK(LRW))
    endif
    ITOL=1
    ITASK=1
    IOPT=1
    ISTATE=1

    IWORK(5:10)=(/4, 100000, 0, 0, 0, 0/)
    RWORK(5:10)=0.0D0
    !RWORK(7) = HMIN
end subroutine

subroutine mylsode(F, x, NEQ, dt, tstart, tstop, atol, rtol)
    implicit none
    integer, intent(in) :: NEQ
    REAL*8, intent(inout) :: x(NEQ)
    real*8, intent(in) :: tstart, tstop, atol, rtol
    real*8, intent(inout) :: dt
    external F, JAC

    CALL DLSODE(F,NEQ,x,tstart,tstop,ITOL,RTOL,ATOL,ITASK,ISTATE, &
                IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
end subroutine

subroutine  JAC()
end subroutine

    
subroutine rk4(F, x, NEQ, dt, nsteps)
    implicit none
    integer, intent(in) :: NEQ
    REAL*8, intent(inout) :: x(NEQ)
    real*8, intent(in) :: dt!, tstart, tstop
    integer, intent(in) :: nsteps
    REAL*8 :: x1(NEQ), xp(NEQ,5)
    integer :: i
  

    do i = 1, nsteps
            call F(NEQ, dt, x, xp(:,1))

            x1 = x + 0.5*dt*xp(:,1)
            call F(NEQ, dt, x1, xp(:,2))

            x1 = x + 0.5*dt*xp(:,2)
            call F(NEQ, dt, x1, xp(:,3))

            x1 = x + dt*xp(:,3)
            call F(NEQ, dt, x1, xp(:,4))

            x = x + (dt/6.0)*(xp(:,1) + 2.0*(xp(:,2) + xp(:,3)) + xp(:,4))
    enddo
end subroutine rk4

subroutine ek(F, x, dt, tstart, tstop,atol, rtol)
    implicit none
    REAL*8, intent(inout) :: x(:)
    real*8, intent(in) :: tstart, tstop, atol, rtol
    real*8, intent(inout) :: dt
    REAL*8 :: x1(size(x)), rmserr, t
    external F

    t = tstart
    x1 = x
    do while (t<tstop)
        call ekstep(F, x1, dt, atol, rtol, rmserr)
        !write (*,*) t, dt, rmserr
        if (rmserr <= 1.0) then
            x=x1
            if (t + dt > tstop) then
                exit
            endif
            t = t + dt
            if (rmserr <=0.1) then
              dt = dt * 1.52799165273419
            endif
        else
            dt = dt/1.90779652734191
            x1 = x
        endif
    enddo
    call ekstep(F, x, tstop-t, atol, rtol, rmserr)
end subroutine


subroutine ekstep(F, x, dt, atol, rtol, rmserr)
    implicit none
    REAL*8, intent(inout) :: x(:)
    real*8, intent(in) :: dt, atol, rtol
    real*8, intent(out) :: rmserr
    integer, parameter :: p=2
    REAL*8 :: x1(size(x)), xp(size(x),p), xe(size(x))
    integer :: NEQ
    external F

    NEQ = size(x)
    call F(NEQ, dt, x, xp(:,1))
    x1 = x + dt*xp(:,1)
    call F(NEQ, dt, x1, xp(:,2))

    x = x + 0.5*dt*(xp(:,1) + xp(:,2))
    xe = 0.5*dt*(xp(:,1) - xp(:,2))

    xe = xe / (abs(rtol*x) + atol)
    rmserr = sqrt( sum(xe**2) / NEQ )
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
    integer :: j, k
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
    real*8 :: t0 = 0.0
    external F

    call F(NEQ, dt, x, xp(:,1))
    do j=2,p
        x1 = x
        do k=1,j-1
                x1 = x1 + a(k,j)*dt*xp(:,k)
        enddo
        call F(NEQ, t0+dt*c(j), x1, xp(:,j))
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
