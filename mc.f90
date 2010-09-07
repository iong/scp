module mc
    use ljmc
    implicit none
    integer, parameter :: nmoves = 4
    real*8, allocatable :: r(:,:), rnew(:,:), rmin(:,:)
    real*8, allocatable :: xstep(:), Z(:), kT(:), beta(:), Cv(:)
    integer, allocatable :: naccepted(:), ntrials(:), tpool(:), stepdim(:)
    integer :: ptinterval
    real*8 :: Tmin, Tmax, U0, Umin
    integer :: ptsynctag = 1, ptswaptag=2, mcblocktag = 3
contains

subroutine pol2cart(pol, cart)
    real*8, intent(in) :: pol(3)
    real*8, intent(out) :: cart(3)
    real*8 :: rxy

    rxy = pol(1) * sin(pol(2))
    cart(1) = rxy * cos(pol(3))
    cart(2) = rxy * sin(pol(3))
    cart(3) = pol(1) * cos(pol(2))
end subroutine

subroutine mc_move_atom(x, rad)
    implicit none
    real*8, intent(inout) :: x(3)
    real*8, intent(in) :: rad
    real*8, dimension(3) :: pol, cart

    pol(1) = gaussran(rad, 0.0d0)
    call random_number(pol(2:3))
    pol(2) = acos(2.0*pol(2) - 1.0)
    pol(3) = 2.0*M_PI*pol(3)
    call pol2cart(pol, cart)
    x = x + cart
end subroutine

subroutine rnd_tuple1(length, nmax, tuple)
    integer, intent(in) :: length, nmax
    integer, intent(out) ::  tuple(length)
    real*8 :: rn
    integer :: i, j
    logical :: unique

    do i=1,length
        do
            unique = .TRUE.
            call random_number(rn)
            tuple(i) = 1 + aint(rn*nmax)
            do j=1,i-1
                if ( tuple(i) == tuple(j) ) then
                    unique = .FALSE.
                    exit
                endif
            enddo
            if (unique) exit
        enddo
    enddo
end subroutine

subroutine rnd_tuple2(length)
    integer, intent(in) :: length
    real*8 :: rn
    integer :: i, r, t, nleft

    nleft = size(tpool)
    do i=1,length
        call random_number(rn)
        r = aint(rn*nleft)
        t=tpool(i+r)
        tpool(i+r) = tpool(i)
        tpool(i) = t
        nleft = nleft - 1
    enddo
end subroutine



subroutine mc_trial(istep)
    use vgw
    implicit none
    integer, intent(in) :: istep
    integer, parameter :: acceptance_trials = 1000
    integer :: i, k, j
    real*8 :: lUmin, Unew, p, rn

    write (*, '(I1$)') me
    rnew(:,:) = r(:,:)
    lUmin = 1d10
    ntrials(istep) = ntrials(istep) + 1
    call rnd_tuple2(stepdim(istep))

    do i=1,stepdim(istep)
        j = tpool(i)
        call mc_move_atom(rnew(:,j), xstep( istep))

        do k=1,3
            if (rnew(k,j) > bl2) then
                rnew(k,j) = rnew(k,j) - bl
            elseif (rnew(k,j) < -bl2) then
                rnew(k,j) = rnew(k,j) + bl
            endif
        enddo
    enddo

    call vgw0(rnew(:,:), Unew, beta(me+1), 0.0d0, y0)
    p = exp( - (Unew - U0) * beta(me+1) )
    call random_number(rn)

    if (p>rn) then
        U0 = Unew
        r(:,:) = rnew(:,:)
        naccepted(istep) = naccepted(istep) + 1
    endif

    if (U0 < Umin) then
        Umin = U0
        rmin(:,:) = r(:,:)
    endif

    if (mod(ntrials(istep), acceptance_trials) == 0) then
        if (naccepted(istep) < 0.3*acceptance_trials) then
            xstep(istep) = xstep(istep) / 1.10779652734191
        elseif (naccepted(istep) > 0.4*acceptance_trials) then
            xstep(istep) = xstep(istep) * 1.12799165273419;
        endif
        naccepted(istep) = 0
    endif
end subroutine

subroutine mc_block(blen, sublen)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: blen, sublen
    integer :: i, j, istep, ierr, req(nprocs), idx, blen2
    real*8 :: rn
    logical :: flag

    blen2 = blen
    if (me /= 0) then
         blen2 = blen*100
    end if
    do i=1,blen2
        flag = .FALSE.
        call MPI_Iprobe(0, mcblocktag, MPI_COMM_WORLD, flag, MPI_STATUS_IGNORE, IERR)
        if (flag) then
            call MPI_Recv(idx, 1, MPI_INTEGER, 0, mcblocktag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            return
        end if

        if (mod(i-1,sublen) == 0) then
            call random_number(rn)
            istep = 1 + aint(nmoves*rn)
        end if

        call mc_trial(istep)
        Z(me+1) = Z(me+1) + exp(-beta(me+1) * U0)
        call pt_swap(i)
    enddo

    if (me==0) then
        do i=1,(nprocs-1)
            call MPI_Isend(1, 1, MPI_INTEGER, i, mcblocktag, MPI_COMM_WORLD, req(i), ierr)
        end do
        do i=1,(nprocs-1)
            call MPI_Waitany(nprocs-1, req(2:nprocs), idx, MPI_STATUS_IGNORE, ierr)
        end do
    end if
end subroutine


subroutine mc_burnin(burnlen)
    integer, intent(in) :: burnlen
    integer ::  j, k

    xstep(1) = 0.1
    do j=1,nmoves
        if (j>1) xstep(j) = xstep(j-1)
        do k=1,burnlen
            call mc_trial(j)
        enddo
    enddo
end subroutine

subroutine pt_swap(nmcsteps)
    use vgw
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nmcsteps
    real*8 :: rbuf(3*Natom), rn, p
    real*8 :: U_rlow(2),U_rhigh(2), betalow, betahigh
    integer :: ilow, ihigh, s(MPI_STATUS_SIZE)
    logical :: flag
    integer :: dest, IERR

    flag = .FALSE.
    call MPI_Iprobe(MPI_ANY_SOURCE, ptsynctag, MPI_COMM_WORLD, flag, s, IERR)
    if (flag) then
            ilow = s(MPI_SOURCE)
            ihigh = me
            call MPI_Recv(rn,1,MPI_REAL8, ilow, ptsynctag, MPI_COMM_WORLD, s, ierr)
    elseif (me<(nprocs-1) .and. mod(nmcsteps,ptinterval)==0) then
            ilow = me
            ihigh = me+1
            call random_number(rn)
            call MPI_Send(rn, 1, MPI_REAL8, ihigh, ptsynctag, MPI_COMM_WORLD, ierr)
    else
        return
    endif

    betahigh = beta(ihigh+1)
    betalow = beta(ilow+1)
    if (me==ilow) then
        U_rlow(1) = U0
        call vgw0(r(:,:), U_rlow(2), betahigh, 0.0d0, y0)
        call MPI_Sendrecv(U_rlow, 2, MPI_REAL8, ihigh, ptswaptag, &
                U_rhigh, 2, MPI_REAL8, ihigh, ptswaptag, MPI_COMM_WORLD, s, IERR)
    else
        call vgw0(r(:,:), U_rhigh(1), betalow, 0.0d0, y0)
        U_rhigh(2) = U0
        call MPI_Sendrecv(U_rhigh, 2, MPI_REAL8, ilow, ptswaptag, &
                U_rlow, 2, MPI_REAL8, ilow, ptswaptag, MPI_COMM_WORLD, s, IERR)
    endif
    
    p = exp(-betalow*(U_rhigh(1) - U_rlow(1))-betahigh*(U_rlow(2) - U_rhigh(2)))
    if (p>rn) then
        write (*,*) 'swapping', ilow,ihigh, 'p=',p
        dest = ilow
        if (me == ilow) dest=ihigh
        call MPI_Sendrecv(reshape(r, (/3*Natom/)), 3*Natom, MPI_REAL8, dest, &
            ptswaptag, &
            rbuf, 3*Natom, MPI_REAL8, dest, ptswaptag, &
            MPI_COMM_WORLD, s, IERR)
        r = reshape(rbuf, (/3, Natom/) )
    endif
end subroutine


end module mc
