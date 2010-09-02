subroutine setup_ljmc()
    use ljmc
    use mc
    implicit none
    integer :: i

    allocate (q0(3,natom), y0(1+21*natom), &
        r(3,natom), rnew(3,natom), &
        rmin(3,natom), ntrials(nmoves), naccepted(nmoves),&
        stepdim(nmoves),xstep(nmoves), tpool(natom), &
        Z(nprocs), kT(nprocs), beta(nprocs), Cv(nprocs))

    kT(1:nprocs) = (/(Tmin + (Tmax-Tmin)/(nprocs-1)*(i-1), i=1,nprocs)/)
    beta(1:nprocs) = 1.0/kT(1:nprocs)

    stepdim=(/(1 + i*(Natom-1)/(nmoves-1),i=0,nmoves-1)/)

    tpool(1:natom) = (/(i,i=1,natom)/)
end subroutine

