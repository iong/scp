module mc
    use ljmc
    implicit none
    integer, parameter :: nmoves = 4, ntau = 10
    real*8, allocatable :: r(:,:), rnew(:,:), rmin(:,:)
    real*8, allocatable :: xstep(:), Z(:), kT(:), beta(:), &
                        Zlocal(:), taugrid(:), U0(:)
    integer, allocatable :: naccepted(:), ntrials(:), tpool(:), stepdim(:)
    integer :: ptinterval
    real*8 :: Tmin, Tmax, Umin
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
    integer, parameter :: acceptance_trials = 1000, ntau = 5
    integer :: i, k, j
    real*8 :: lUmin, Unew(ntau), taugrid(ntau), p, rn

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
    integer :: i, j, istep, ierr, req(nprocs), NMC2, nmclast, nmcmaster, s(MPI_STATUS_SIZE)
    real*8 :: rn
    logical :: flag

    NMC2 = NMC
    if (me /= 0) then
         NMC2 = NMC*100
    end if
    Zlocal = 0.0
    nmclast = 0
    do i=1,NMC2
         if (mod(i-1,sublen) == 0) then
            call random_number(rn)
            istep = 1 + aint(nmoves*rn)
        end if

        call mc_trial(istep)
        Zlocal = Zlocal + exp(-taugrid * (U0 - U0(ntau)))

        call MPI_Iprobe(MPI_ANY_SOURCE, ptsynctag, MPI_COMM_WORLD, flag, s, IERR)
        if (flag) then
            call pt_swap(s(MPI_SOURCE), me)
        end if

        call MPI_Iprobe(0, mcblocktag, MPI_COMM_WORLD, flag, MPI_STATUS_IGNORE, IERR)
        if (flag) then
            call MPI_Recv(nmcmaster, 1, MPI_INTEGER, 0, mcblocktag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            call mc_dump_state(i, nmclast, nmcmaster)
            nmclast = i
            Zlocal = 0.0
            if (nmcmaster >= NMC) then
                return
            end if
        end if

        if (me==0 .and. mod(i,mcblen) == 0) then
            do j=1,(nprocs-1)
                call MPI_Isend(i, 1, MPI_INTEGER, j, mcblocktag, MPI_COMM_WORLD, req(j+1), ierr)
            end do
            call MPI_Waitall(nprocs-1, req(2:nprocs), MPI_STATUS_IGNORE, ierr)
            call mc_dump_state(i, nmclast, i)
            nmclast = i
            Zlocal = 0.0
            if (nmcmaster >= NMC) then
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
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nmcnow, nmclast, nmcmaster
    real*8 :: Cvlocal(ntau)
    integer :: ierr, i
    character(256) :: fname, label
    character(20) :: pestr, nowstr, masterstr
 

    write(label, "('U0 =',F14.7)") U0(ntau)
    write(fname, "('dump/pe',I3,'/dump_',I10,'_',I10,'.xyz')") me, nmcmaster, nmcnow
    call replace_char(fname, ' ', '0')
    call dump_xyz(r, fname, label)

    Zlocal = Zlocal / real(nmcnow-nmclast, 8)
    call heat_capacity(ntau, Zlocal, taugrid, Cvlocal)

    write(fname, "('dump/pe',I3,'/Z_',I10,'_',I10,'.dat')") me, nmcmaster, nmcnow
    call replace_char(fname, ' ', '0')
    open(30, file=fname,STATUS='REPLACE')
    write (30,'(4(ES16.8," "))') (1.0/taugrid(i), Cvlocal(i), Z(i), beta(i),i=ntau,1,-1)
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
        call vgw0s(r(:,:), U_rlow(2), betahigh, 0.0d0, y0)
        call MPI_Sendrecv(U_rlow, 2, MPI_REAL8, ihigh, ptswaptag, &
                U_rhigh, 2, MPI_REAL8, ihigh, ptswaptag, MPI_COMM_WORLD, s, IERR)
    else
        call vgw0s(r(:,:), U_rhigh(1), betalow, 0.0d0, y0)
        U_rhigh(2) = U0(ntau)
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
