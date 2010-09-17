module mc
    use ljmc
    implicit none
    integer, parameter :: nmoves = 4, ntau = 10
    real*8, allocatable :: r(:,:), rnew(:,:), rmin(:,:)
    real*8, allocatable :: xstep(:), kT(:), beta(:), &
                        Z(:), taugrid(:), U0(:)
    integer, allocatable :: naccepted(:), ntrials(:), tpool(:), stepdim(:)
    integer :: ptinterval, dumpinterval
    real*8 :: Tmin, Tmax, Umin
    integer :: ptsynctag = 1, ptswaptag=2, dumptag = 3
    integer :: ptswaps = 0, ptattempts = 0
contains

subroutine mc_move_atom(x, rad)
    use utils
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
    use utils
    implicit none
    integer, intent(in) :: istep
    integer, parameter :: acceptance_trials = 1000
    integer :: i, k, j
    real*8 :: lUmin, Unew(ntau), p, rn

    !write (*, '(I1$)') me
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

    call vgw0(rnew(:,:), Unew, taugrid, 0.0d0, y0)
    p = exp( - (Unew(ntau) - U0(ntau)) * taugrid(ntau) )
    call random_number(rn)

    if (p>rn) then
        U0 = Unew
        r(:,:) = rnew(:,:)
        naccepted(istep) = naccepted(istep) + 1
        !if (me/=nprocs-1 .and. abs(U0(ntau) - U0(1)) > 100) then
        !    write (*,*) 'On PE', me
        !    write (*,"((2F14.7))") (taugrid(i), U0(i),i=1,ntau)
        !endif
    endif

    if (U0(ntau) < Umin) then
        Umin = U0(ntau)
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


subroutine mc_run_loop(NMC, mcblen, sublen)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: NMC, mcblen, sublen
    integer :: i, j, istep, ierr, req(nprocs), NMC2, nmclast, nmcmaster
    integer :: s(MPI_STATUS_SIZE), statuses(MPI_STATUS_SIZE,nprocs)
    real*8 :: rn
    logical :: flag

    NMC2 = NMC
    if (me /= 0) then
         NMC2 = NMC*100
    end if
    Z = 0.0
    nmclast = 0
    do i=1,NMC2
         if (mod(i-1,sublen) == 0) then
            call random_number(rn)
            istep = 1 + aint(nmoves*rn)
        end if

        call mc_trial(istep)
        Z = Z + exp(-taugrid * (U0 - U0(ntau)))

        call MPI_Iprobe(MPI_ANY_SOURCE, ptsynctag, MPI_COMM_WORLD, flag, s, IERR)
        if (flag) then
            call pt_swap(s(MPI_SOURCE), me)
        end if

        call MPI_Iprobe(0, dumptag, MPI_COMM_WORLD, flag, MPI_STATUS_IGNORE, IERR)
        if (flag) then
            call MPI_Recv(nmcmaster, 1, MPI_INTEGER, 0, dumptag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            call mc_dump_state(i, nmclast, nmcmaster)
            if (mod(nmcmaster,mcblen) == 0) then 
                nmclast = i
                Z = 0.0
            end if
            if (nmcmaster >= NMC) then
                return
            end if
        end if

        if (me==0 .and. mod(i,dumpinterval) == 0) then
            do j=1,(nprocs-1)
                call MPI_Isend(i, 1, MPI_INTEGER, j, dumptag, MPI_COMM_WORLD, req(j), ierr)
            end do
            call MPI_Waitall(nprocs-1, req, statuses, ierr)
            call mc_dump_state(i, nmclast, i)
            if (mod(i,mcblen) == 0) then 
                nmclast = i
                Z = 0.0
            end if
            if (i >= NMC) then
                return
            end if
        endif

        if (me<(nprocs-1) .and. mod(i,ptinterval)==0) then
            call pt_swap(me, me+1)
        end if
    enddo
end subroutine

subroutine mc_dump_state(nmcnow, nmclast, nmcmaster)
    use xyz
    use utils
    use ljmc
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nmcnow, nmclast, nmcmaster
    real*8 :: Cv(ntau)
    integer :: ierr, i
    character(256) :: fname, label

    write(label, "('U0 =',F14.7)") U0(ntau)
    write(fname, "('dump/pe',I3,'/dump_',I10,'.xyz')") me, nmcmaster
    call replace_char(fname, ' ', '0')
    call dump_xyz(r, fname, label)

    call heat_capacity2(Z(1:ntau), taugrid(1:ntau), Cv(1:ntau))

    write(fname, "('dump/pe',I3,'/Z_',I10,'.dat')") me, nmcmaster
    call replace_char(fname, ' ', '0')
    open(30, file=fname)
    write (30,'(4(ES16.8," "))') (1.0/taugrid(i), Cv(i), Z(i), taugrid(i),i=ntau,1,-1)
    close(30)

    write(fname, "('dump/pe',I3,'/state_',I10,'.dat')") me, nmcmaster
    call replace_char(fname, ' ', '0')
    open(30, file=fname)
    write (30, "('nmcnow =', I10)") nmcnow
    write (30, "('xstep =', (F8.5))") xstep
    close(30)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
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

subroutine pt_swap(ilow, ihigh)
    use vgw
    use utils
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ilow, ihigh
    real*8 :: rbuf(3*Natom), rn, p
    real*8 :: U_rlow(2),U_rhigh(2), betalow, betahigh
    integer :: s(MPI_STATUS_SIZE)
    integer :: dest, IERR

    if (me==ilow) then
        call random_number(rn)
        call MPI_Send(rn, 1, MPI_REAL8, ihigh, ptsynctag, MPI_COMM_WORLD, ierr)
        ptattempts = ptattempts + 1
    else if (me==ihigh) then
        call MPI_Recv(rn,1,MPI_REAL8, ilow, ptsynctag, MPI_COMM_WORLD, s, ierr)
    else
        write (*,*) 'me is neither ilow or ihigh'
        stop
    endif

    betahigh = beta(ihigh+1)
    betalow = beta(ilow+1)

    if (me==ilow) then
        U_rlow(1) = U0(ntau)
        call vgw0s(r, U_rlow(2), betahigh, 0.0d0, y0)
        !call vgw0(r, U_rlow, (/betahigh, betalow/), 0.0d0, y0)
        !call fliplr(U_rlow)
        call MPI_Sendrecv(U_rlow, 2, MPI_REAL8, ihigh, ptswaptag, &
                U_rhigh, 2, MPI_REAL8, ihigh, ptswaptag, MPI_COMM_WORLD, s, IERR)
    else
        call vgw0s(r, U_rhigh(1), betalow, 0.0d0, y0)
        U_rhigh(2) = U0(ntau)
        !call vgw0(r, U_rhigh, (/betahigh, betalow/), 0.0d0, y0)
        !call fliplr(U_rhigh)
        call MPI_Sendrecv(U_rhigh, 2, MPI_REAL8, ilow, ptswaptag, &
                U_rlow, 2, MPI_REAL8, ilow, ptswaptag, MPI_COMM_WORLD, s, IERR)
    endif
    
    p = exp(-betalow*(U_rhigh(1) - U_rlow(1))-betahigh*(U_rlow(2) - U_rhigh(2)))
    !if (me==ilow) write (*,*) me, 'p =', p
    !write (*,*) betalow, betahigh, U_rhigh, U_rlow, U0(ntau)
    if (p>rn) then
        dest = ilow
        if (me == ilow) then
            ptswaps = ptswaps + 1
            write (*,'("swapping ",I2,"-",I2,5X,"p =",F8.5)') ilow, ihigh, real(ptswaps)/real(ptattempts)
            dest=ihigh
        end if
        call MPI_Sendrecv(reshape(r, (/3*Natom/)), 3*Natom, MPI_REAL8, dest, &
            ptswaptag, &
            rbuf, 3*Natom, MPI_REAL8, dest, ptswaptag, &
            MPI_COMM_WORLD, s, IERR)
        r = reshape(rbuf, (/3, Natom/) )
        call vgw0(r, U0, taugrid, 0.0d0, y0)
    endif
end subroutine


end module mc
