program lj
    use xyz
    use utils
   ! use propagation
    implicit none
    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer, parameter :: eout = 30, cvvout = 31, track_width = 3
    integer :: Natom, Nbath, ntracks, tracksep, seglen, ringpos = 1
    real*8 ::tstart, tstop, dt, bl, bl2, rho, kT
    real*8, dimension(:), allocatable :: Qbath, xi, vxi
    real*8, dimension(:,:), allocatable :: r0,  r, p, f, rshift
    real*8, dimension(:,:,:), allocatable :: rt0, pt0, track, r0shift, rprev, pprev
    integer, allocatable :: trackstart(:), timestamp(:)
    character(256) :: stem

    integer :: iostat
    integer :: ndt, ne, nequil
    real*8 ::Ekin ,Epot,tequil, tsep, tlen, pcm(3), Q1nhc
    character(LEN=256) :: cfgfile, coords
    integer :: i, j, k, n, ixyz

    real*8 :: mass
    namelist /gmdcfg/mass,kT,rho,coords,tstart,tstop,dt,Nbath,Q1nhc,ne,tequil,tlen,tsep
    
    interface
        subroutine nose_hoover_chain(p, Ekin, kT, xi, vxi, Q, dt, ne)
            real*8, intent(inout) :: p(:,:), Ekin, xi(:), vxi(:)
            real*8, intent(in) :: kT, Q(:), dt
            integer, intent(in) :: ne
        end subroutine nose_hoover_chain
    end interface

    cfgfile='lj.in'

    open(20, file=cfgfile)
    read(20, NML=gmdcfg)
    close(20)

    if (command_argument_count() >0) then
       call get_command_argument(1, coords)
    end if

    call load_xyz(coords, r0)
    Natom = size(r0, 2)
    
    ndt = (tstop - tstart) / dt
    nequil = tequil / dt
    tracksep = tsep / dt
    seglen = tlen / dt
    ntracks = tlen / tsep

    allocate (r(3,natom), p(3,natom), f(3,natom), rshift(3,Natom), &
            rt0(3,natom,ntracks), pt0(3,natom,ntracks), r0shift(3,natom,ntracks), &
            track(track_width,seglen,ntracks), trackstart(ntracks), &
            rprev(3,Natom,seglen), pprev(3,Natom,seglen), &
            timestamp(seglen))
    call seed_rng()

    bl = (Natom/rho)**(1.0/3.0)
    bl2=bl/2

    mass = mass*0.020614788876D0

    if (Nbath>=2) then
        allocate(Qbath(Nbath),  xi(Nbath), vxi(Nbath))
        Qbath = kT*Q1nhc
        Qbath(1) = 3*Natom*Qbath(1)
        do i=1,Nbath
            vxi(i) = gaussran(sqrt(kT/Qbath(i)), 0.0d0)
        end do
        xi = 0.0d0
        vxi = 0.0d0
    end if

    do i=1,Natom
        do j=1,3
            p(j,i) = gaussran(sqrt(kT*mass), 0.0d0)
        end do
    end do
    pcm =  sum(p, 2)/Natom
    do i=1,Natom
        p(:,i) = p(:,i) - pcm
    end do

    write (*,*) shape(r0shift)
    r = r0
    pt0 = 0.0d0
    rt0 = 0.0d0
    r0shift = 0.0d0
    rshift = 0.0d0
    track = 0.0d0
    trackstart = (/ (tracksep*i + seglen, i=0,ntracks-1) /)
    write (*,*) 'y'
    Ekin = 0.5*sum(p**2)/mass
    if (Nbath > 0) then
        call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, 0.0d0, ne)
    end if
    call lj_f(r, f, Epot)
    if (Nbath > 0) then
        call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, 0.0d0, ne)
    end if


    ixyz = index(coords, '.xyz', .TRUE.)
    stem = 'dump/lj'//coords(1:ixyz-1)
    open(eout,file=trim(stem)//'_energy.dat')

    write(eout,'(6F18.7)') 0.0d0,Ekin, Epot,Ekin+Epot
    do i=1,ndt
        if (Nbath > 0) then
            call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, dt, ne)
        end if
        p = p + 0.5*dt*f
        r = r + (dt/mass)*p

        call lj_f(r, f, epot)
        p = p + 0.5*dt*f

        Ekin = 0.5*sum(p**2)/mass
        if (Nbath > 0) then
            call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, dt, ne)
        end if
        
        if (i>nequil) then
            call correlations(i-nequil)
        end if

        do j=1,Natom
            do k=1,3
                if (abs(r(k, j)) > bl2) then
                    rshift(k, j) = rshift(k, j) + sign(bl, r(k, j))
                    r(k, j) = r(k, j) - sign(bl, r(k, j))
                end if
            end do
        end do
        write(eout,'(6F18.7)') i*t0fs,Ekin, Epot,Ekin+Epot
    end do
contains

subroutine lj_f(r, f, U)
    real*8, intent(in) :: r(:,:)
    real*8, intent(out) :: f(:,:), U
    integer :: i, j, Natom
    real*8 :: fij(3), dr(3), y6, r2
    real*8, parameter :: eps0 = 32.35d0, sigma0 = 3.045d0


    Natom = size(r, 2)
    f = 0.0d0
    U = 0.0d0
    do i=1,Natom-1
        do j=i+1,Natom
            dr = r(:,j) - r(:,i)
            where (abs(dr) > bl2)
                dr = dr - sign(bl, dr)
            end where
            r2 = sum(dr**2)
            y6 = (sigma0*sigma0/r2)**3
            fij = 4.0d0*eps0*(12.0d0*y6**2-6.0d0*y6)*dr/r2
            U = U + 4.0d0*eps0*(y6*y6 - y6)

            f(:,j) = f(:,j) + fij
            f(:,i) = f(:,i) - fij
        end do
    end do
end subroutine

subroutine correlations(ndt)
    use utils
    implicit none
    integer, intent(in) :: ndt
    integer :: i,j, trackpos

    pprev(:,:,ringpos) = p
    rprev(:,:,ringpos) = r + rshift
    timestamp(ringpos) = ndt

    do i=1,ntracks
        if (ndt == trackstart(i)) then
            call init_track(i, ndt)
        end if
    end do
    
    do i=1,ntracks-1
        if (trackstart(i) > ndt) cycle
        if (ndt - trackstart(i)  == seglen) then
            track(:,:,i+1) = track(:,:,i+1) + track(:,:,i)
            call init_track(i, ndt)
        end if
    end do

    if (ndt - trackstart(ntracks) == seglen) then
        track(:,:,ntracks) = track(:,:,ntracks) / (2*Natom * ntracks)
        call dump_track(track(:,:,ntracks), ndt/seglen - 1)
        call init_track(ntracks, ndt)
    end if
    
    do i=1,ntracks
        if (trackstart(i) > ndt) cycle
        trackpos = ndt - trackstart(i) + 1
        call update_track(i, trackpos, r+rshift, p)
    end do
    ringpos  = mod(ringpos, seglen) + 1
end subroutine

subroutine init_track(trackno, ndt)
    implicit none
    integer, intent(in) :: trackno, ndt
    integer :: j

    rt0(:,:,trackno) = r
    pt0(:,:,trackno) = p
    r0shift(:,:,trackno) = rshift

    track(:,:,trackno) = 0.0d0
    trackstart(trackno) = ndt

    do j=1,seglen
        call update_track(trackno, ndt - timestamp(j) + 1, &
                rprev(:,:,j), pprev(:,:,j))
    end do
    track(2, :, trackno) = -track(2, :, trackno)
end subroutine init_track

subroutine update_track(trackno, trackpos, r, p)
    integer, intent(in) :: trackno, trackpos
    real*8, dimension(3,Natom), intent(in) :: r, p

    track(1, trackpos, trackno) =  track(1, trackpos, trackno) + sum(pt0(:,:,trackno)*p)/mass**2
    track(2, trackpos, trackno) =  track(2, trackpos, trackno) + sum((rt0(:,:,trackno) + r0shift(:,:,trackno))*p)/mass
    track(3, trackpos, trackno) =  track(3, trackpos, trackno) + sum((rt0(:,:,trackno) + r0shift(:,:,trackno))*r)
end subroutine update_track


subroutine dump_track(tr, trackno)
    use utils
    implicit none
    real*8, intent(in) :: tr (:,:)
    integer, intent(in) :: trackno
    integer :: i
    character(4) :: cdump

    call int2strz(trackno, 4, cdump)
    open(cvvout, file=trim(stem)//'_Cvv_'//cdump//'.dat')
    write(cvvout,'(4F18.7)') (dt*(i-1)*t0fs, tr(1:track_width,i), i=1,seglen)
    close(cvvout)
end subroutine


end program lj
! vim:et
