module kubo_mod
    use utils
    use mainvars
contains
subroutine kubo(q0, qv0, beta, nsteps, xk, vk)

    implicit none
    real*8, intent(in) :: q0(:,:), qv0(:,:), beta
    integer, intent(in) :: nsteps
    real*8, intent(out) :: xk(:,:), vk(:,:)
    real*8, allocatable, target :: yl(:,:)
    real*8, pointer:: ql(:), qr(:), Qnkl(:), Qnkr(:)
    real*8, dimension(3*fm%Natom, 3*fm%Natom) :: gl, gr, glrinv, glr
    real*8, dimension(3*fm%Natom) :: gqv, xksum, vksum, vnkl, vnkr
    real*8 :: tau, dtau, w
    integer :: i, qGQnkLength, qnkpos, info
 
    if (mod(nsteps,2) /= 0) then
        write (*,*) 'nsteps must be even!!'
        stop
    end if

    qGQnkLength = 3*fm%Natom + fm%NG + fm%NQnk
    qnkpos = 3*fm%Natom + fm%NG

    allocate (yl(qGQnkLength , 0:nsteps) )

    tau = 0.0d0
    dtau = beta / nsteps

    call fm%Ueff(q0, 0d0)
    yl(:,0) = fm%y ( 1 : qGQnkLength )

    do i=1,nsteps
        tau = tau + dtau
        call fm%propagate(tau)

        yl(:,i) = fm%y(1 : qGQnkLength)
    end do

    xksum = yl(1:3*fm%Natom, nsteps/2)
    Qnkl => yl(qnkpos +1 : qnkpos + fm%NQnk, nsteps/2)
    call dsymv('U', 3*fm%Natom, 1d0, Qnkl, 3*fm%Natom, qv0, 1, 0d0, vksum, 1)

    do i=0,nsteps/2-1

        ql => yl(1:3*fm%Natom, i)
        call fm_get_g(fm%Natom, yl(:,i), gl)
        Qnkl => yl(qnkpos + 1 : qnkpos + fm%NQnk, i)

        qr => yl(1:3*fm%Natom, nsteps - i)
        call fm_get_g(fm%Natom, yl(:,nsteps - i), gr)
        Qnkr => yl(qnkpos + 1 : qnkpos + fm%NQnk, nsteps - i)

        ! trapeze integration
        if (i==0) then
            w = 1.0d0
        else
            w = 2.0d0
        end if

        glrinv = gl + gr
        glr = glrinv
        call dpotrf('U', 3 * fm%Natom, glrinv, 3*fm%Natom, info)
        call dpotri('U', 3 * fm%Natom, glrinv, 3*fm%Natom, info)

        call dsymv('U', 3*fm%Natom, 1d0, gr, 3*fm%Natom, ql, 1, 0d0, gqv, 1)
        if (maxval(abs(gqv)) > 1e20) then
            call break_in_debugger()
        end if
        call dsymv('U', 3*fm%Natom, 1d0, gl, 3*fm%Natom, qr, 1, 1d0, gqv, 1)
        if (maxval(abs(gqv)) > 1e20) then
            call break_in_debugger()
        end if
        vnkl = xksum
        call dsymv('U', 3*fm%Natom, w, glrinv, 3*fm%Natom, gqv, 1, 1d0, xksum, 1)
        xksum = xksum + vnkr*w
        if (maxval(abs(xksum)) > 1e20) then
            call break_in_debugger()
        end if

        call dsymv('U', 3*fm%Natom, 1d0, Qnkl, 3*fm%Natom, qv0, 1, 0d0, vnkl, 1)
        call dsymv('U', 3*fm%Natom, 1d0, Qnkr, 3*fm%Natom, qv0, 1, 0d0, vnkr, 1)

        call dsymv('U', 3*fm%Natom, 1d0, gr, 3*fm%Natom, vnkl, 1, 0d0, gqv, 1)
        if (maxval(abs(gqv)) > 1e20) then
            call break_in_debugger()
        end if
        call dsymv('U', 3*fm%Natom, 1d0, gl, 3*fm%Natom, vnkr, 1, 1d0, gqv, 1)
        if (maxval(abs(gqv)) > 1e20) then
            call break_in_debugger()
        end if

        call dsymv('U', 3*fm%Natom, w, glrinv, 3*fm%Natom, gqv, 1, 1d0, vksum, 1)
        if (maxval(abs(vksum)) > 1e20) then
            call break_in_debugger()
        end if
    end do
        if (maxval(abs(xksum)) > 1e20) then
            call break_in_debugger()
        end if

    xk = reshape(xksum / nsteps, (/ 3, fm%Natom /) )
    vk = reshape(vksum / nsteps, (/ 3, fm%Natom /) )

    deallocate (yl)
end subroutine

subroutine break_in_debugger()
    write (*,*) 'Ooooooooo!'
end subroutine
end module kubo_mod
