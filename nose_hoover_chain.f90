subroutine nose_hoover_chain(p, Ekin, kT, xi, vxi, Q, dt, ne)
!
! Molecular Physics 87, 1117 (1996)
!
    real*8, intent(inout) :: p(:,:), Ekin, xi(:), vxi(:)
    real*8, intent(in) :: kT, Q(:), dt
    integer, intent(in) :: ne
    integer, parameter :: nys = 3
    real*8, dimension(3),parameter :: w = (/1.35120719195966d0, -1.70241438391932d0, 1.35120719195966d0/)
    real*8 :: dts, GM, Gk, G1, scale
    integer :: i, j, k, M, Nf

    M = size(xi)
    scale = 1.0d0
    Nf = size(p)
    do i=1,ne
        do j=1,nys
            dts = dt * w(j) / ne
            GM = (Q(M-1)*vxi(M-1)**2 - kT)/Q(M)
            vxi(M) = vxi(M) + 0.025d0 * dts * GM
            
            do k=M-1,2,-1
                vxi(k) = vxi(k) * exp( -0.125d0 * dts * vxi(k+1) )
                Gk = (Q(k-1) * vxi(k-1)**2 - kT) / Q(k)
                vxi(k) = vxi(k) + 0.25d0 * dts * Gk
                vxi(k) = vxi(k) * exp( -0.125d0 * dts * vxi(k+1) )
            end do
                
            vxi(1) = vxi(1) * exp( -0.125d0 * dts * vxi(2) )
            G1 = (2.0d0*Ekin - Nf*kT) / Q(1)
            vxi(1) = vxi(1) + 0.25d0 * dts * G1
            vxi(1) = vxi(1) * exp( -0.125d0 * dts * vxi(2) )

            scale = scale * exp( -0.5d0 * dts * vxi(1) )
            Ekin = Ekin * exp( - dts * vxi(1) )

            xi = xi + dts*vxi

            vxi(1) = vxi(1) * exp( -0.125d0 * dts * vxi(2) )
            G1 = (2.0d0*Ekin - Nf*kT) / Q(1)
            vxi(1) = vxi(1) + 0.25d0 * dts * G1
            vxi(1) = vxi(1) * exp( -0.125d0 * dts * vxi(2) )

            do k=2,M-1
                vxi(k) = vxi(k) * exp( -0.125d0 * dts * vxi(k+1) )
                Gk = (Q(k-1) * vxi(k-1)**2 - kT) / Q(k)
                vxi(k) = vxi(k) + 0.25d0 * dts * Gk
                vxi(k) = vxi(k) * exp( -0.125d0 * dts * vxi(k+1) )
            end do

            dts = dt * w(j) / ne
            GM = (Q(M-1)*vxi(M-1)**2 - kT)/Q(M)
            vxi(M) = vxi(M) + 0.025d0 * dts * GM
        end do
    end do

    p = p * scale
end subroutine
