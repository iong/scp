
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
    use ljmc
    implicit none
    real*8, intent(inout) :: x(3)
    real*8, intent(in) :: rad
    real*8, dimension(3) :: pol, cart

    pol(1) = gaussran(rad, 0.0)
    call random_number(pol(2:3))
    pol(2) = acos(2.0*pol(2) - 1.0)
    pol(3) = 2.0*M_PI*pol(3)
    call pol2cart(pol, cart)
    x = x + cart
end subroutine

subroutine mc_1by1(mcburn, irep)
    use ljmc
    implicit none
    integer, intent(in) :: mcburn, irep
    integer, parameter :: acceptance_trials = 1000
    integer :: i, j, k
    real*8 :: lUmin, Unew, p, rn,bl2
    rnew(:,:,irep) = r(:,:,irep)
    call vgw0(rnew(:,:,irep), U0(irep), beta(irep), 0.0, y0)
    lUmin = 1e10
    bl2 = bl/2
    do i=1,mcburn
        ntrials(irep) = ntrials(irep) + 1
        call random_number(rn)
        j = 2 + aint((Natom-1)*rn)
        call mc_move_atom(rnew(:,j,irep), rmove(irep))


        do k=1,3
            if (rnew(k,j,irep) > bl2) then
                rnew(k,j,irep) = rnew(k,j,irep) - bl
            elseif (rnew(k,j,irep) < -bl2) then
                rnew(k,j,irep) = rnew(k,j,irep) + bl
            endif
        enddo

        call vgw0(rnew(:,:,irep), Unew, beta(irep), 0.0, y0)
        p = exp( - (Unew - U0(irep)) * beta(irep) )
        call random_number(rn)

        if (p>rn) then
            U0(irep) = Unew
            r(:,:,irep) = rnew(:,:,irep)
            naccepted(irep) = naccepted(irep) + 1
        endif

        if (U0(irep) < Umin(irep)) then
            Umin(irep) = U0(irep)
            rmin(:,:,irep) = r(:,:,irep)
        endif

        if (mod(ntrials(irep), acceptance_trials) == 0) then
            if (naccepted(irep) < 0.3*acceptance_trials) then
                rmove(irep) = rmove(irep) / 1.10779652734191
            elseif (naccepted(irep) > 0.4*acceptance_trials) then
                rmove(irep) = rmove(irep) * 1.12799165273419;
            endif
            naccepted(irep) = 0
        endif
    enddo
end subroutine

subroutine pt_swap()
    use ljmc
    implicit none
    real*8 :: rtmp(3,Natom), rn, p
    real*8 :: U_betai_ri,U_betai_rj,U_betaj_ri,U_betaj_rj
    integer :: i, j
    call random_number(rn)
    i = 1 + aint( rn * (nstreams-1) )
    j = i+1

    U_betai_ri = U0(i);
    U_betaj_rj = U0(j);
    call vgw0(r(:,:,i), U_betaj_ri, beta(j), 0.0, y0)
    call vgw0(r(:,:,j), U_betai_rj, beta(i), 0.0, y0)

    p = exp(-beta(i)*(U_betai_rj - U_betai_ri)-beta(j)*(U_betaj_ri - U_betaj_rj))
    call random_number(rn)
    if (p>rn) then
        rtmp(:,:) = r(:,:,i)
        r(:,:,i) = r(:,:,j)
        r(:,:,j) = rtmp(:,:)

    endif
    return
end subroutine
