subroutine correlations(ndt)
    use spine
    use vgw
    use utils
    implicit none
    integer, intent(in) :: ndt
    integer :: i,j, trackpos
    real * 8 :: vs(3,Natom)

    do j=1,Natom
        v(:,j) = matmul(invMeff(:,:,j), p(:,j))
        vs(:,j) = matmul(matsqrt(invMeff(:,:,j)), p(:,j))/sqrt(mass)
    end do

    do i=1,ntracks
        if (ndt == trackstart(i)) then
            call init_track(i, ndt, vs)
        end if
    end do
    
    do i=1,ntracks-1
        if (trackstart(i) > ndt) cycle
        if (ndt - trackstart(i)  == seglen) then
            track(:,:,i+1) = track(:,:,i+1) + track(:,:,i)
            call init_track(i, ndt, vs)
        end if
    end do

    if (ndt - trackstart(ntracks) == seglen) then
        track(:,:,ntracks) = track(:,:,ntracks) / (Natom * ntracks)
        call dump_track(track(:,:,ntracks), ndt/seglen)
        call init_track(ntracks, ndt, vs)
    end if

    do i=1,ntracks
        if (trackstart(i) > ndt) cycle
        trackpos = ndt - trackstart(i) + 1

        track(1, trackpos, i) =  track(1, trackpos, i) + sum(v0tau(:,:,i)*v)
        track(2, trackpos, i) =  track(2, trackpos, i) + sum(v0s(:,:,i)*vs)
        track(3, trackpos, i) =  track(3, trackpos, i) + sum(p0(:,:,i)*v)
        track(4, trackpos, i) =  track(4, trackpos, i) + sum(vkubo(:,:,i)*v)
    end do
end subroutine

subroutine init_track(trackno, ndt, vs)
    use vgw
    use spine
    implicit none
    integer, intent(in) :: trackno, ndt
    real*8, intent(in) :: vs(3,Natom)
    real*8 :: xk(3,Natom)
    integer :: j

    call unpack_Qnk(y, Qnk)
    do j=1,Natom
        v0tau(:,j,trackno) = matmul(Qnk(:,:,j), v(:,j))
    end do
    v0s(:,:,trackno) = vs
    p0(:,:,trackno) = p/mass
    call kubo(r, v, 1.0d0/kT, 200, xk, vkubo(:,:,trackno))
    write(31,*)  p0(:,:,trackno)
    write(32,*)  vkubo(:,:,trackno)

    track(:,:,trackno) = 0.0d0
    trackstart(trackno) = ndt
end subroutine init_track

subroutine dump_track(tr, trackno)
    use spine
    use utils
    implicit none
    real*8, intent(in) :: tr (:,:)
    integer, intent(in) :: trackno
    integer :: i
    character(4) :: cdump

    call int2strz(trackno, 4, cdump)
    open(cvvout, file=trim(stem)//'_Cvv_'//cdump//'.dat')
    write(cvvout,'(5F18.7)') (dt*(i-1)*t0fs, tr(1:4,i), i=1,seglen)
    close(cvvout)
end subroutine


subroutine mean_sq_disp(dxsq)
    use spine
    real*8, intent(out) :: dxsq

    dxsq = sum((r - r0equil - rshift)**2)/Natom
end subroutine
