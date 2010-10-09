module spine
    implicit none
    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer :: Natom, Nbath
    real*8 :: rcmin, tstart, tstop, dtout, dt, bl, bl2,rho, kT
    real*8, dimension(:), allocatable :: qp, y, ekin, epot, etot, Cvv, Qbath
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
    real*8 :: Ueff, v3(3), p3(3), Ekin
    integer :: i, Nq, Nqp


    Nq = 3*Natom
    Nqp = 6*Natom
    call vgw1(reshape(x(1:Nq), (/3, Natom/)), Ueff, 1.0/kT, 0.0d0, y, Meff, invMeff)
    lastepot = -2.0*kT*y(1)

    Ekin = 0
    do i=1,Natom
        p3 = x(Nq+3*i-2:Nq+3*i)
        v3 = matmul(invMeff(:,:,i), p3)
        Ekin = Ekin + dot_product(v3, p3)
        xp(3*i-2:3*i) = v3
    end do

    xp(1+Nq:Nqp) = 2.0*kT*y(2+18*Natom:1+21*Natom)

    if (Nbath <= 0) then
        return
    end if

    xp(1+Nq:Nqp) = xp(1+Nq:Nqp) - x(1+Nq:Nqp)*x(Nqp+1)/Qbath(1)

    xp(Nqp+1) = Ekin - Nq*kT - x(Nqp+1)*x(Nqp+2)/Qbath(2)
    do i=2,Nbath-1
        xp(Nqp+i) = x(Nqp+i-1)**2/Qbath(i-1) - kT - x(Nqp+i)*x(Nqp+i+1)/Qbath(i+1)
    end do
    xp(Nqp+Nbath) = x(Nqp+Nbath-1)**2/Qbath(Nbath-1) - kT
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

    etot(i) = ekin(i) + epot(i)
end subroutine

end module spine
