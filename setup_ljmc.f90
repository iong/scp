subroutine setup_ljmc()
    use ljmc
    use mc
    implicit none
    integer :: i

    allocate (q0(3,natom), y0(1+21*natom), &
        r(3,natom,nstreams), rnew(3,natom,nstreams), &
        rmin(3,natom,nstreams), U0(nstreams), &
        Umin(nstreams), ntrials(nstreams,nmoves), naccepted(nstreams,nmoves),&
        stepdim(nmoves),xstep(nstreams,nmoves), tpool(natom), &
        Z(nstreams), kT(nstreams), beta(nstreams), Cv(nstreams))

    kT(1:nstreams) = (/(Tmin + (Tmax-Tmin)/(nstreams-1)*(i-1), i=1,nstreams)/)
    beta(1:nstreams) = 1.0/kT(1:nstreams)

    stepdim=(/(1 + i*(Natom-1)/(nmoves-1),i=0,nmoves-1)/)

    tpool(1:natom) = (/(i,i=1,natom)/)
end subroutine

