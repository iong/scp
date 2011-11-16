module correlations
    use mainvars
    use utils
    implicit none
contains
subroutine update_track(ndt)
    implicit none
    integer, intent(in) :: ndt
    integer :: j
    real * 8 :: vs(3,Natom), rs(3,Natom), sqrtMeff(3,3), rt(3,Natom), &
        d0k(3), dt(3)

    call dsymv('U', 3*Natom, 1d0, invMeff, 3*Natom, p, 1, 0d0, v, 1)
    call dsymv('U', 3*Natom, 1d0, sqrtInvMeff, 3*Natom, p, 1, 0d0, vs, 1)
    call dsymv('U', 3*Natom, 1d0, sqrtMeff, 3*Natom, r+rshift, 1, &
        0d0, rs, 1)

    do j=1,Natom
        vs(:,j) = vs(:,j)/sqrt(fm%mass(j))
        rs(:,j) = rs(:,j)/sqrt(fm%mass(j))
    end do

    rt =  r + rshift
    !write(*,*) ndt, track(1, ndt), sum( v0tau * v)/Natom

    d0k = r0k(:,1) - r0k(:,2)
    dt = rt(:,1) - rt(:,2)

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
    track(11, ndt) = track(11, ndt) +  sum( d0k * dt)
    track(12, ndt) = track(12, ndt) +  sum( r0s *   rs)
    track(13, ndt) = track(13, ndt) +  sum( v0s *   rs)
    track(14, ndt) = track(14, ndt) +  sum( r0s *   vs)
end subroutine

subroutine init_track()
    use kubo_mod
    implicit none
    integer :: j, yskip
    double precision, pointer :: pQnk(:)

    call dgemv('N', 3*Natom, 3*Natom, 1d0, invMeff, 3*Natom, p, 1, 0d0, v, 1)
    call dgemv('N', 3*Natom, 3*Natom, 1d0, sqrtInvMeff, 3*Natom, p, 1, 0d0, v0s, 1)
    call dgemv('N', 3*Natom, 3*Natom, 1d0, sqrtMeff, 3*Natom, r+rshift, 1, &
        0d0, r0s, 1)

    do j=1,Natom
        v0(:,j) = p(:,j)/fm%mass(j)
        v0s(:,j) = v0s(:,j)/sqrt(fm%mass(j))
        r0s(:,j) = r0s(:,j)/sqrt(fm%mass(j))
    end do

    q0tau = fm%get_q() + rshift

    yskip = 3*fm%Natom + fm%NG 
    pQnk => fm%y(yskip + 1 : yskip + fm%NQnk)

    call dgemv('N', 3*Natom, 3*Natom, 1d0, pQnk, 3*Natom, v, 1, 0d0, v0tau, 1)

    call kubo(r, v, 1.0d0/kT, 200,  r0k, vkubo)
    r0shift = rshift
    v0k =  0.0d0*(r0k - rkold) / dt

    r0k = r0k + rshift

    track = 0.0d0
end subroutine init_track


subroutine dump_track(tr, trackno)
    use utils
    implicit none
    real*8, intent(in) :: tr (:,:)
    integer, intent(in) :: trackno
    integer :: i
    character(4) :: cdump

    call int2strz(trackno, 4, cdump)
    write(*,*) trim(stem)//'_Cvv_'//cdump//'.dat'
    open(cvvout, file=trim(stem)//'_Cvv_'//cdump//'.dat')
    write(cvvout,'(15F18.7)') (dt*(i-1), tr(:,i), i=1,ndt)
    close(cvvout)

    open(mapout, file='dump_map')
    write(mapout, '("kT"/"CvvSym"/"EasyKubo1"/"EasyKubo2"/"CvvKubo"/"q0tau_v"/"r0k_v"/"q0tau_r")')
    write(mapout, '("r0k_r"/"v0tau_r"/"dr0k_r"/"vkubo_r"/"r0s_rs"/"v0s_rs"/"r0s_vs")')
    close(mapout)

end subroutine
end module correlations
