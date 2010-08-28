subroutine heat_capacity(Z, kT, Cv)
    use ljmc
    implicit none
    real*8, intent(in), dimension(nstreams) :: Z, kT
    real*8, intent(out), dimension(nstreams) :: Cv
    integer :: j
    real*8, dimension(nstreams) :: dT, Tmid, dlogZ

    dT = kT(2:nstreams) - kT(1:nstreams-1)
    Tmid = 0.5*(kT(2:nstreams) + kT(1:nstreams-1))
    dlogZ = log(Z(2:nstreams) / Z(1:nstreams-1))

    dlogZdT = Tmid**2 * dlogZ/dT
    do j=2,nstreams-1
        Cv(i) = (dlogZ(j) - dlogZ(j-1)) *2.0/ (T(j+1)-T(j-1))
    enddo

    Cv(1) = Cv(2)
    Cv(nstreams) = Cv(nstreams-1)
end subroutine
