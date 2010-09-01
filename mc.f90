module mc
    use ljmc
    implicit none
    integer, parameter :: nmoves = 4
    real*8, allocatable :: r(:,:,:), rnew(:,:,:), rmin(:,:,:), U0(:), Umin(:)
    real*8, allocatable :: xstep(:,:), Z(:), kT(:), beta(:), Cv(:)
    integer, allocatable :: naccepted(:,:), ntrials(:,:), tpool(:), stepdim(:)
    integer :: nstreams, ptinterval
    real*8 :: Tmin, Tmax
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



subroutine mc_trial(irep, istep)
    use vgw
    implicit none
    integer, intent(in) :: irep, istep
    integer, parameter :: acceptance_trials = 100
    integer :: i, k, j
    real*8 :: lUmin, Unew, p, rn

    rnew(:,:,irep) = r(:,:,irep)
    lUmin = 1d10
    ntrials(irep,istep) = ntrials(irep,istep) + 1
    call rnd_tuple2(stepdim(istep))

    do i=1,stepdim(istep)
        j = tpool(i)
        call mc_move_atom(rnew(:,j,irep), xstep(irep, istep))

        do k=1,3
            if (rnew(k,j,irep) > bl2) then
                rnew(k,j,irep) = rnew(k,j,irep) - bl
            elseif (rnew(k,j,irep) < -bl2) then
                rnew(k,j,irep) = rnew(k,j,irep) + bl
            endif
        enddo
    enddo

    call vgw0(rnew(:,:,irep), Unew, beta(irep), 0.0d0, y0)
    p = exp( - (Unew - U0(irep)) * beta(irep) )
    call random_number(rn)

    if (p>rn) then
        U0(irep) = Unew
        r(:,:,irep) = rnew(:,:,irep)
        naccepted(irep,istep) = naccepted(irep,istep) + 1
    endif

    if (U0(irep) < Umin(irep)) then
        Umin(irep) = U0(irep)
        rmin(:,:,irep) = r(:,:,irep)
    endif

    if (mod(ntrials(irep,istep), acceptance_trials) == 0) then
        if (naccepted(irep,istep) < 0.3*acceptance_trials) then
            xstep(irep,istep) = xstep(irep,istep) / 1.10779652734191
        elseif (naccepted(irep,istep) > 0.4*acceptance_trials) then
            xstep(irep,istep) = xstep(irep,istep) * 1.12799165273419;
        endif
        naccepted(irep,istep) = 0
    endif
end subroutine

subroutine pt_swap()
    use vgw
    implicit none
    real*8 :: rtmp(3,Natom), rn, p
    real*8 :: U_betai_ri,U_betai_rj,U_betaj_ri,U_betaj_rj
    integer :: i, j
    call random_number(rn)
    i = 1 + aint( rn * (nstreams-1) )
    j = i+1

    U_betai_ri = U0(i);
    U_betaj_rj = U0(j);
    call vgw0(r(:,:,i), U_betaj_ri, beta(j), 0.0d0, y0)
    call vgw0(r(:,:,j), U_betai_rj, beta(i), 0.0d0, y0)

    p = exp(-beta(i)*(U_betai_rj - U_betai_ri)-beta(j)*(U_betaj_ri - U_betaj_rj))
    call random_number(rn)
    if (p>rn) then
        write (*,*) 'swapping', i,j, 'p=',p
        rtmp(:,:) = r(:,:,i)
        r(:,:,i) = r(:,:,j)
        r(:,:,j) = rtmp(:,:)
    endif
    return
end subroutine

subroutine mc_block(blen, sublen)
    integer, intent(in) :: blen, sublen
    integer :: irep, nsub, i, j, k, istep
    real*8 :: rn

    nsub = blen/sublen
    do i=1,nsub
        call random_number(rn)
        istep = 1 + aint(nmoves*rn)
        do j=1,sublen,ptinterval
            do irep=1,nstreams
                do k=1,ptinterval
                    call mc_trial(irep, istep)
                    Z(irep) = Z(irep) + exp(-beta(irep) * U0(irep))
                enddo
            enddo
            call pt_swap()
        enddo
    enddo
end subroutine


subroutine mc_burnin(burnlen)
    integer, intent(in) :: burnlen
    integer :: irep, j, k

    xstep(:,1) = 0.1
    do irep=1,nstreams
        do j=1,nmoves
            if (j>1) xstep(irep, j) = xstep(irep, j-1)
            do k=1,burnlen
                call mc_trial(irep, j)
            enddo
        enddo
    enddo
end subroutine

end module mc
