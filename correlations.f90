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

    vprev(:,:,ringpos) = v
    vsprev(:,:,ringpos) = vs
    rprev(:,:,ringpos) = r + rshift
    timestamp(ringpos) = ndt

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
        track(:,:,ntracks) = track(:,:,ntracks) / (2*Natom * ntracks)
        call dump_track(track(:,:,ntracks), ndt/seglen - 1)
        call init_track(ntracks, ndt, vs)
    end if
    
    do i=1,ntracks
        if (trackstart(i) > ndt) cycle
        trackpos = ndt - trackstart(i) + 1
        call update_track(i, trackpos, r+rshift, v, vs)
    end do
    ringpos  = mod(ringpos, seglen) + 1
end subroutine

subroutine init_track(trackno, ndt, vs)
    use vgw
    use spine
    implicit none
    integer, intent(in) :: trackno, ndt
    real*8, intent(in) :: vs(3,Natom)
    real*8 :: Ueff
    integer :: j


    call vgw1(r, Ueff, 1.0d0/kT, 0.0d0, y)
    call unpack_Qnk(y, Qnk)
    do j=1,Natom
        v0tau(:,j,trackno) = matmul(Qnk(:,:,j), v(:,j))
    end do
    v0s(:,:,trackno) = vs
    p0(:,:,trackno) = p/mass


    call unpack_q(y, q0tau(:,:,trackno))
    call kubo(r, v, 1.0d0/kT, 200,  r0k(:,:,trackno), vkubo(:,:,trackno))
    r0shift(:,:,trackno) = rshift

    track(:,:,trackno) = 0.0d0
    trackstart(trackno) = ndt

    do j=1,seglen
        call update_track(trackno, ndt - timestamp(j) + 1, &
                rprev(:,:,j), vprev(:,:,j), vsprev(:,:,j))

    end do
    track((/ 7, 8 /), :, trackno) = -track((/ 7, 8 /), :, trackno)
end subroutine init_track

subroutine update_track(trackno, trackpos, r, v, vs)
    use spine
    integer, intent(in) :: trackno, trackpos
    real*8, dimension(3,Natom), intent(in) :: r, v, vs

    track(1, trackpos, trackno) =  track(1, trackpos, trackno) + sum(v0tau(:,:,trackno)*v)
    track(2, trackpos, trackno) =  track(2, trackpos, trackno) + sum(v0s(:,:,trackno)*vs)
    track(3, trackpos, trackno) =  track(3, trackpos, trackno) + sum(p0(:,:,trackno)*v)
    track(4, trackpos, trackno) =  track(4, trackpos, trackno) + sum(vkubo(:,:,trackno)*v)
    track(5, trackpos, trackno) =  track(5, trackpos, trackno) + sum(q0tau(:,:,trackno)*(r - r0shift(:,:,trackno)))
    track(6, trackpos, trackno) =  track(6, trackpos, trackno) + sum(r0k(:,:,trackno)*(r - r0shift(:,:,trackno)))
    track(7, trackpos, trackno) =  track(7, trackpos, trackno) + sum(q0tau(:,:,trackno)*v)
    track(8, trackpos, trackno) =  track(8, trackpos, trackno) + sum(r0k(:,:,trackno)*v)
end subroutine update_track


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
    write(cvvout,'(9F18.7)') (dt*(i-1)*t0fs, tr(1:track_width,i), i=1,seglen)
    close(cvvout)
end subroutine
