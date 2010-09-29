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
       
include 'potential_energy.f90'
include 'vgw0.f90'
include 'vgw1.f90'
include 'rhss0.f90'
include 'rhss1.f90'

end module vgw

