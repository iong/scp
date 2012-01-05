program OH1D
    use correlations
    use mainvars
    use utils
    use vgwfm_mod
    use xyz
    implicit none
    integer :: iostat

    integer ::np=100, nx=51, ip, ix, ipb
    double precision :: beta, xeq, dx, w0, Ekin, Z, Epot, Epotref, pcm(3), &
          mw_tanh_hwbeta

    character(LEN=256) :: argin, fname
    integer :: i, j, k, ixyz

    double precision, parameter :: A_au = 1.8897, K_au = 3.16681534282952d-06, &
        kOH = 0.49536d0

    Natom = 2


    if (command_argument_count() >0) then
       call get_command_argument(1, argin)
       read(argin,*) kT
    end if
    if (command_argument_count() >1) then
       call get_command_argument(2, argin)
       read (argin, *) nx
    end if

    kT = kT * K_au
    beta = 1d0/kT

    xeq = 0.9824 * A_au !equilibrium of the quantum OH bond
    mw_tanh_hwbeta = sqrt(1823.0*kOH)*tanh(0.5d0*sqrt(kOH/1823.0)/kT)
    dx = 3.5d0 /sqrt(mw_tanh_hwbeta)

    dt = 2.0*M_PI*sqrt(1823.0/kOH)/32.0


    tstart = 0d0
    tstop = 32*5*dt
    ndt = (tstop - tstart) / dt

    allocate (fm, rkold(3, Natom), &
                Meff(3*natom, 3*Natom), invMeff(3*natom, 3*natom), &
                sqrtMeff(3*natom, 3*Natom), sqrtInvMeff(3*natom, 3*natom), &
                r(3,natom), r0(3,Natom), p(3,natom), p0(3,natom), v(3,natom), &
                f(3, Natom), &
                q0tau(3,natom), rshift(3,natom),r0k(3,natom), &
                v0k(3,natom), r0shift(3,natom), &
                v0tau(3,natom), v0s(3,natom), r0s(3,natom), &
                v0(3,natom), vkubo(3,natom), &
                track(1:track_width,0:ndt), &
                trackaccum(1:track_width,0:ndt))

    call seed_rng()

    call fm%init(2, 'OH')
    call fm%set_mm(.TRUE.)

    trackaccum = 0.0d0

    r0 = 0d0
    Z = 0d0
    do ix = 1,nx
        r0(1,2) = (2d0*real(ix-1)/real(nx-1) - 1d0)*dx + xeq
        write (*,*) 'ix = ', ix, nx
        do ip=1,2*np
            r = r0
            call fm%Ueff(r, 1.0/kT, Epot)
            call fm%get_Meff(Meff)

            w0 = exp(-beta*Epot) * r0(1,2)**2
            if (mod(ip,2) == 1) then
                call initial_momenta(kT, Meff, p0)
                p = p0
            else
                p = -p0
            end if

            call verletstep(0.0d0, Epot)
            call init_track()
            call update_track(0)
        
           
            Ekin = 0.5*sum(p * v)
            if (debug .and. ip==1) then
                write(eout,'(6F18.7)') 0.0d0,Ekin, Epot,Ekin+Epot
            end if
            do i=1,ndt
                call verletstep(dt, Epot)
                call update_track(i)
                Ekin = 0.5*sum(p * v)
                if (debug .and. ip==1) then
                    write(eout,'(6F18.7)') dt*i,Ekin, Epot,Ekin+Epot
                end if
            end do
            if (debug .and. ip==1) then
                write (eout, '(//)')
            end if
            trackaccum = trackaccum + track*w0
            Z = Z + w0
        end do
    end do
    trackaccum = trackaccum / (2*Natom*Z) 
    write (*,*) trackaccum(1,0)
    stem='OH'
    call dump_track(trackaccum, 0)
    close(eout)

contains


subroutine initial_momenta(kT, Meff, p)
    implicit none
    real*8, intent(in) :: kT, Meff(:,:)
    real*8, intent(out) :: p(:, :)
    real*8 :: U(size(Meff, 1),size(Meff, 1))
    double precision, dimension (size(Meff, 1)) :: W, FV1, FV2, pr
    integer :: i, j, ierr, N
    
    N = size(Meff, 1)
    call RS(N,N,Meff,W,1,U,FV1,FV2,ierr)
    do j=1,N
        pr(j) = gaussran(sqrt(W(j)*kT), 0.0d0)
    end do
    call dgemv('N', N, N, 1d0, U, N, pr, 1, 0d0, p, 1)
end subroutine


subroutine verletstep(dt, Epot)
    implicit none
    real*8, intent(in) :: dt
    real*8, intent(out) :: Epot
    integer :: i, k

    p = p + 0.5*dt*f
    call dsymv('U', 3*Natom, dt, invMeff, 3*Natom, p, 1, 1d0, r, 1)

    call fm%Ueff(r, 1.0/kT, Epot, f)
    call fm%get_Meff(Meff, invMeff, sqrtMeff, sqrtInvMeff)

    p = p + 0.5*dt*f
end subroutine

end program OH1D
! vim:et
