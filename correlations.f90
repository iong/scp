subroutine correlations(ndt)
    use spine
    use vgw
    use utils
    implicit none
    integer, intent(in) :: ndt
    integer :: i,j, trackpos
    real * 8 :: vs(3,Natom), rs(3,Natom), sqrtMeff(3,3), vkwaste(3,Natom)

    
    do j=1,Natom
        sqrtMeff = matsqrt(Meff(:,:,j))
        v(:,j) = matmul(invMeff(:,:,j), p(:,j))
        vs(:,j) = matmul(matsqrt(invMeff(:,:,j)), p(:,j))/sqrt(mass)
        rs(:,j) = matmul(sqrtMeff, r(:,j) +  rshift(:,j))/sqrt(mass)
    end do

    vprev(:,:,ringpos) = v
    vsprev(:,:,ringpos) = vs
    rprev(:,:,ringpos) = r + rshift
    rsprev(:,:,ringpos) = rs
    
    timestamp(ringpos) = ndt

    do i=1,ntracks
        if (ndt == trackstart(i)) then
            call init_track(i, ndt, rs, vs)
        end if
        if (ndt - trackstart(i)  == seglen - 1) then
            call kubo(r, v, 1.0d0/kT, 200,  rkold, vkwaste)
        end if
    end do

    

    do i=1,ntracks-1
        if (trackstart(i) > ndt) cycle
        if (ndt - trackstart(i)  == seglen) then
            track(:,:,i+1) = track(:,:,i+1) + track(:,:,i)
            call init_track(i, ndt, rs, vs)
        end if
    end do

    if (ndt - trackstart(ntracks) == seglen) then
        trackaccum = trackaccum + track(:,:,ntracks) / (2*Natom * ntracks)
        call dump_track(trackaccum / (ndt/seglen), ndt/seglen - 1)
        call init_track(ntracks, ndt, rs, vs)
    end if
    
    do i=1,ntracks
        if (trackstart(i) > ndt) cycle
        trackpos = ndt - trackstart(i) + 1
        call update_track(i, trackpos, r+rshift, rs, v, vs)
    end do
    ringpos  = mod(ringpos, seglen) + 1
end subroutine

subroutine init_track(trackno, ndt, rs, vs)
    use vgw
    use spine
    implicit none
    integer, intent(in) :: trackno, ndt
    real*8, intent(in) :: rs(3,Natom), vs(3,Natom)
    real*8 :: Ueff
    integer :: j


    call unpack_Qnk(y, Qnk)
    do j=1,Natom
        v0tau(:,j,trackno) = matmul(Qnk(:,:,j), v(:,j))
    end do
    v0s(:,:,trackno) = vs
    p0(:,:,trackno) = p/mass

    r0s(:,:,trackno) = rs
    call unpack_q(y, q0tau(:,:,trackno))
    call kubo(r, v, 1.0d0/kT, 200,  r0k(:,:,trackno), vkubo(:,:,trackno))
    r0shift(:,:,trackno) = rshift
    v0k(:,:,trackno) =  (r0k(:,:,trackno) - rkold) / dt

    q0tau(:,:,trackno) = q0tau(:,:,trackno) + rshift
    r0k(:,:,trackno) = r0k(:,:,trackno) + rshift

   

    track(:,:,trackno) = 0.0d0
    trackstart(trackno) = ndt

    do j=1,seglen
        call update_track(trackno, ndt - timestamp(j) + 1, &
                rprev(:,:,j), rsprev(:,:,j), &
                vprev(:,:,j), vsprev(:,:,j))

    end do
    track((/ 5,6,9,10,11,13,14 /), :, trackno) = -track((/ 5,6,9,10,11,13,14 /), :, trackno)
end subroutine init_track

subroutine update_track(trackno, trackpos, rt, rst, vt, vst)
    use spine
    integer, intent(in) :: trackno, trackpos
    real*8, dimension(3,Natom), intent(in) :: rt, rst, vt, vst
    real*8 :: trackslice(track_width)

    trackslice(1) = sum(  v0tau(:,:,trackno) * vt)
    trackslice(2) = sum(  v0s(:,:,trackno) *   vst)
    trackslice(3) = sum(  p0(:,:,trackno) *    vt)
    trackslice(4) = sum(  vkubo(:,:,trackno) * vt)
    trackslice(5) = sum(  q0tau(:,:,trackno) * vt) 
    trackslice(6) = sum(  r0k(:,:,trackno) *   vt)
    trackslice(7) = sum(  q0tau(:,:,trackno) * rt)
    trackslice(8) = sum(  r0k(:,:,trackno) *   rt)
    trackslice(9) = sum(  v0tau(:,:,trackno) * rt)
    trackslice(10) = sum( v0k(:,:,trackno) *   rt)
    trackslice(11) = sum( vkubo(:,:,trackno) * rt)
    trackslice(12) = sum( r0s(:,:,trackno) *   rst)
    trackslice(13) = sum( v0s(:,:,trackno) *   rst)
    trackslice(14) = sum( r0s(:,:,trackno) *   vst)

    track(:,trackpos, trackno) = track(:,trackpos, trackno) + trackslice
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
    !open(cvvout, file=trim(stem)//'_Cvv_'//cdump//'.dat')
    open(cvvout, file=trim(stem)//'_Cvv.dat')
    write(cvvout,'(15F18.7)') (dt*(i-1)*t0fs, tr(1:track_width,i), i=1,seglen)
    close(cvvout)
end subroutine
