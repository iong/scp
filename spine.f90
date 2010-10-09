module gmd_mod
    implicit none
    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer :: Natom
    real*8 :: rcmin, tstart, tstop, dtout, dt, bl, bl2,rho, kT
    real*8, dimension(:), allocatable :: qp, y, ekin, epot, etot, Cvv
    real*8, dimension(:,:), allocatable :: r0, p0, vtau0, v0
    real*8, dimension(:,:,:), allocatable :: Meff0, invMeff0, Qnk0, Meff, invMeff
    real*8 :: lastepot
contains
subroutine initial_momenta(kT, Meff, p)
    use utils
    implicit none
    real*8, intent(in) :: kT, Meff(:,:,:)
    real*8, intent(out) :: p(:, :)
    real*8 :: U(3,3), W(3), FV1(3), FV2(3), p3(3)
    integer :: i, j, ierr
    
    do i=1,size(Meff, 3)
        call RS(3,3,Meff(:,:,i),W,1,U,FV1,FV2,ierr)
        do j=1,3
            p3(j) = gaussran(sqrt(W(j)*kT), 0.0d0)
        end do
        p(:,i) = matmul(U, p3)
    end do
end subroutine


subroutine RHS_Veff(NEQ, dt, x, xp)
    use vgw
    implicit none
    integer, intent(in) :: NEQ
    real*8, intent(in) :: dt, x(NEQ)
    real*8, intent(out) :: xp(NEQ)
    real*8 :: Ueff
    integer :: i, Natom

    Natom = NEQ/6

    call vgw1(reshape(x(1:NEQ/2), (/3, Natom/)), Ueff, 1.0/kT, 0.0d0, y, Meff, invMeff)

    do i=1,Natom
        xp(1+3*(i-1):3*i) = matmul(invMeff(:,:,i), x(1+3*(Natom+i-1):3*(Natom+i)))
    end do

    xp(1+3*Natom:6*Natom) = 2.0*kT*y(2+18*Natom:1+21*Natom)
    lastepot = -2.0*kT*y(1)
end subroutine 

subroutine update_TCF(i)
    implicit none
    integer, intent(in) :: i
    real*8 :: v3(3), p3(3)
    integer :: j

    ekin(i) = 0
    do j=1,Natom
        p3 = qp(1+3*(Natom+j-1):3*(Natom+j))
        v3 = matmul(invMeff(:,:,j), p3)
        Cvv(i) = Cvv(i) + dot_product(vtau0(:,j), v3)
        ekin(i) = ekin(i) + dot_product(p3, v3)
        epot(i) = lastepot
    end do
    ekin(i) = 0.5*ekin(i)

end subroutine

end module gmd_mod
