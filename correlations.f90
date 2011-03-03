module correlations
    use vgw
    use spine
    use utils
    implicit none
contains
subroutine update_track(ndt)
    implicit none
    integer, intent(in) :: ndt
    integer :: i,j, trackpos
    real * 8 :: vs(3,Natom), rs(3,Natom), sqrtMeff(3,3), rt(3,Natom)
    real*8 :: trackslice(track_width)

    do j=1,Natom
        sqrtMeff = matsqrt(Meff(:,:,j))
        v(:,j) = matmul(invMeff(:,:,j), p(:,j))
        vs(:,j) = matmul(matsqrt(invMeff(:,:,j)), p(:,j))/sqrt(mass)
        rs(:,j) = matmul(sqrtMeff, r(:,j) +  rshift(:,j))/sqrt(mass)
    end do

    rt =  r + rshift
    !write(*,*) ndt, track(1, ndt), sum( v0tau * v)/Natom

    track(1, ndt)  = track(1, ndt)  +  sum( v0tau * v)
    track(2, ndt)  = track(2, ndt)  +  sum( v0s *   vs)
    track(3, ndt)  = track(3, ndt)  +  sum( v0 *    v)
    track(4, ndt)  = track(4, ndt)  +  sum( vkubo * v)
    track(5, ndt)  = track(5, ndt)  +  sum( q0tau * v) 
    track(6, ndt)  = track(6, ndt)  +  sum( r0k *   v)
    track(7, ndt)  = track(7, ndt)  +  sum( q0tau * rt)
    track(8, ndt)  = track(8, ndt)  +  sum( r0k *   rt)
    track(9, ndt)  = track(9, ndt)  +  sum( v0tau * rt)
    track(10, ndt) = track(10, ndt) +  sum( v0k *   rt)
    track(11, ndt) = track(11, ndt) +  sum( vkubo * rt)
    track(12, ndt) = track(12, ndt) +  sum( r0s *   rs)
    track(13, ndt) = track(13, ndt) +  sum( v0s *   rs)
    track(14, ndt) = track(14, ndt) +  sum( r0s *   vs)
end subroutine

subroutine init_track()
    implicit none
    real*8 :: Ueff
    integer :: j
    real * 8 :: sqrtMeff(3,3)

    do j=1,Natom
        sqrtMeff = matsqrt(Meff(:,:,j))
        v(:,j) = matmul(invMeff(:,:,j), p(:,j))
        v0s(:,j) = matmul(matsqrt(invMeff(:,:,j)), p(:,j))/sqrt(mass)
        r0s(:,j) = matmul(sqrtMeff, r(:,j) +  rshift(:,j))/sqrt(mass)
    end do

    call unpack_Qnk(y, Qnk)
    do j=1,Natom
        v0tau(:,j) = matmul(Qnk(:,:,j), v(:,j))
    end do
    v0 = p/mass

    call unpack_q(y, q0tau)
    call kubo(r, v, 1.0d0/kT, 200,  r0k, vkubo)
    r0shift = rshift
    v0k =  (r0k - rkold) / dt

    q0tau = q0tau + rshift
    r0k = r0k + rshift

    track = 0.0d0
end subroutine init_track


subroutine dump_track(tr, trackno)
    use spine
    use utils
    implicit none
    real*8, intent(in) :: tr (:,:)
    integer, intent(in) :: trackno
    integer :: i
    character(4) :: cdump

     write(*,*) 0, tr(1, 0)
    call int2strz(trackno, 4, cdump)
    open(cvvout, file=trim(stem)//'_Cvv_'//cdump//'.dat')
    write(cvvout,'(15F18.7)') (dt*t0fs*(i-1), tr(:,i), i=1,ndt)
    close(cvvout)
end subroutine
end module correlations
