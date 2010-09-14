subroutine setup_ljmc()
    use ljmc
    use mc
    use utils
    implicit none
    integer :: i
    character(3) :: pestr
    character(40) :: cmdstr

    allocate (q0(3,natom), y0(1+21*natom), &
        r(3,natom), rnew(3,natom), &
        rmin(3,natom), ntrials(nmoves), naccepted(nmoves),&
        stepdim(nmoves),xstep(nmoves), tpool(natom), &
        Zlocal(ntau), U0(ntau), taugrid(ntau), kT(nprocs), beta(nprocs))

    call linspace(Tmin, Tmax, nprocs, kT)
    beta(1:nprocs) = 1.0/kT(1:nprocs)
    if (me == nprocs -1) then
        call linspace(0.5*beta(me+1), beta(me+1), ntau, taugrid)
    else
        call linspace(beta(me+2), beta(me+1), ntau-1, taugrid(2:ntau))
        taugrid(1) = 2.0*taugrid(2) - taugrid(3)
    end if

    stepdim=(/(1 + i*(Natom-1)/(nmoves-1),i=0,nmoves-1)/)

    tpool(1:natom) = (/(i,i=1,natom)/)

    call int2strz(me, 3, pestr)
    write(cmdstr, "('mkdir -p dump/pe',A3)") pestr
    call system(cmdstr)
end subroutine

