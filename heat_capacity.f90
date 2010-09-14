subroutine heat_capacity(N, Z, kT, Cv)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in), dimension(N) :: Z, kT
    real*8, intent(out), dimension(N) :: Cv
    integer :: j
    real*8, dimension(N) :: dT, Tmid, dlogZ, dlogZdT

    dT(1:N-1) = kT(2:N) - kT(1:(N-1))
    Tmid(1:N-1) = 0.5*(kT(2:N) + kT(1:N-1))
    dlogZ(1:N-1) = log(Z(2:N) / Z(1:N-1))

    dlogZdT = Tmid**2 * dlogZ/dT
    do j=2,N-1
        Cv(j) = (dlogZdT(j) - dlogZdT(j-1)) *2.0/ (kT(j+1)-kT(j-1))
    enddo

    Cv(1) = Cv(2)
    Cv(N) = Cv(N-1)
end subroutine

subroutine heat_capacity2(N, Z, beta, Cv)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in) :: Z(:), beta(:)
    real*8, intent(out) :: Cv(:)
    real*8, dimension(N) :: dbeta, dlogZ, dlogZ_dbeta, d2logZ_dbeta2

    dbeta(1:N-1) = beta(2:N) - beta(1:N-1)
    dlogZ(1:N-1) = log(Z(2:N) / Z(1:N-1))

    dlogZ_dbeta(1:N-1) = dlogZ(1:N-1)/dbeta(1:N-1)
    d2logZ_dbeta2 = 2.0 * (dlogZ_dbeta(2:N-1) - dlogZ_dbeta(1:N-2))/(beta(3:N) - beta(1:N-2))

    Cv(2:N-1) = beta(2:N-1)**2*d2logZ_dbeta2
    Cv(1) = Cv(2)
    Cv(N) = Cv(N-1)
end subroutine

