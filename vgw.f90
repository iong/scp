module vgw
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    real*8, private :: BL
    logical, allocatable :: QRC(:)
    real*8 :: invmass, RC, ATOL, RTOL, TAUMIN, mass
    integer :: N_atom

contains

subroutine vgwinit(natom, boxlen)
    implicit none
    integer, intent(in) :: natom
    real*8, intent(in) :: boxlen
    N_atom = natom
    BL = boxlen
    invmass=1.0D0/mass

    allocate(QRC(natom*natom))
end subroutine

subroutine unpack_Qnk(y, Qnk)
    real*8, intent(in) :: y(:)
    real*8, intent(out) :: Qnk(:,:,:)
    integer :: i, j, k, CNT, N

    N = size(Qnk, 3)

    CNT = 2+9*N
    DO I=1,N
        DO J=1,3
            DO K=1,3
                QNK(J,K,I)=Y(CNT)
                CNT=CNT+1
            ENDDO
        ENDDO
    ENDDO
end subroutine

subroutine unpack_q(y, q)
    real*8, intent(in) :: y(:)
    real*8, intent(out) :: q(:,:)
    integer :: Natom

    Natom = size(q,2)
    q = reshape(y(2:1+3*Natom), (/3, Natom/) )
end subroutine

subroutine unpack_g(y, G)
    real*8, intent(in) :: y(:)
    real*8, intent(out) :: G(:,:,:)
    integer :: i, j, k, cnt, N

    N = size(G, 3)
    cnt = 1+3*N+1
    do i=1,N
        do j=1,3
            do k=j,3
                G(k, j, i) = y(cnt)
                G(j, k, i) = y(cnt)
                cnt = cnt + 1
            enddo
        enddo
    enddo
end subroutine

subroutine unpack_f(y, kT, f)
    real*8, intent(in) :: y(:), kT
    real*8, intent(out) :: f(:,:)
    integer :: Natom

    Natom = size(f,2)
    f = 2.0*kT*reshape(y(2+18*Natom:1+21*Natom), (/3, Natom/) )
end subroutine

subroutine init_gaussians(q0, tau, y)
    use propagation
    implicit none
    REAL*8, intent(in) :: Q0(3,N_atom), tau
    REAL*8, intent(inout) :: Y(:)
    real*8 :: ULJ, G0(6), E3(9)
    integer :: NY

    NY = size(y)

    call Upot_tau0(Q0,ULJ)

    Y(1) = -tau*ULJ
    Y(2:1+3*N_atom)=reshape(Q0, (/3*N_atom/))

    G0=tau*invmass*(/1.0, 0.0, 0.0, 1.0, 0.0, 1.0/)
    Y(2+3*N_atom:1+9*N_atom) = reshape(spread(G0, 2, N_atom), (/ 6*N_atom /) )

    E3 = (/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/)
    if (NY == 1 + 21*N_atom) then
        Y(2+9*N_atom:1+18*N_atom) = reshape( &
                spread(E3, 2, N_atom), (/ 9*N_atom /) )

        Y(2+18*N_atom:1+21*N_atom) = reshape(Ux_tau0(Q0), (/3*N_atom/) ) * tau
    end if
end subroutine
       
include 'potential_energy.f90'
include 'vgw0.f90'
include 'vgw1.f90'
include 'rhss0.f90'
include 'rhss1.f90'

end module vgw

