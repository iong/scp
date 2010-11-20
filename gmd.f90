program gmd
    use spine
    use vgw
    use xyz
    use utils
    use propagation
    implicit none
    integer :: iostat
    integer :: ndt, ne, nequil
    real*8 :: Ekin ,Epot,tequil, tsep, tlen, pcm(3), Q1nhc, Ueff

    character(LEN=256) :: cfgfile, coords
    integer :: i, j, k, ixyz
    namelist /gmdcfg/Natom,mass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,kT,rho, &
            coords,tstart,tstop,dt,Nbath,Q1nhc,ne,tequil,tlen,tsep


    cfgfile='pH2.in'
    taumin = 1.0e-6
    rtol=1e-4
    atol=1.0e-4
    rc=10.0
    LJA=0.0d0
    LJC=0.0d0

    open(20, file=cfgfile)
    read(20, NML=gmdcfg)
    close(20)

    if (command_argument_count() >0) then
       call get_command_argument(1, coords)
    end if


    ndt = (tstop - tstart) / dt
    nequil = tequil / dt
    tracksep = tsep / dt
    seglen = tlen / dt
    ntracks = tlen / tsep

    allocate (y(1+21*natom), r0(3,natom), rkold(3, Natom), &
                Qnk(3,3,natom), Meff(3,3,natom), invMeff(3,3,natom), &
                r(3,natom), p(3,natom), v(3,natom), f(3, Natom), &
                q0tau(3,natom,ntracks), rshift(3,natom),r0k(3,natom,ntracks), &
                v0k(3,natom,ntracks), r0shift(3,natom,ntracks), &
                v0tau(3,natom,ntracks), v0s(3,natom,ntracks), r0s(3,natom,ntracks), &
                p0(3,natom,ntracks), vkubo(3,natom,ntracks), &
                track(track_width,seglen,ntracks), trackstart(ntracks), &
                rprev(3,Natom,seglen),  rsprev(3,Natom,seglen), &
                vprev(3,Natom,seglen), vsprev(3,Natom,seglen), &
                timestamp(seglen))
    call load_xyz(r0, coords)
    call seed_rng()

    bl = (Natom/rho)**(1.0/3.0)
    bl2=bl/2

    mass = mass*0.020614788876D0
    call vgwinit(natom, bl)

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

    p = 0.0d0
    call vgw1(r0, Ueff, 1.0/kT, 0.0d0, y, Meff, invMeff)
    !call kubo(r0, p,  1.0/kT, 200, r, p)
    !stop
    call initial_momenta(kT, Meff, p)
    pcm = sum(p, 2)/Natom
    do i=1,Natom
        p(:,i) = p(:,i) - pcm
    end do

    r = r0
    v0tau = 0.0d0
    vkubo = 0.0d0
    q0tau = 0.0d0
    r0shift = 0.0d0
    rshift = 0.0d0
    track = 0.0d0
    trackstart = (/ (tracksep*i + seglen, i=0,ntracks-1) /)

    Ekin = kinetic_energy()
    if (Nbath > 0) then
        call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, 0.0d0, ne)
    end if
    call verletstep(0.0d0, Epot)
    if (Nbath > 0) then
        call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, 0.0d0, ne)
    end if


    ixyz = index(coords, '.xyz', .TRUE.)
    stem = 'dump/'//coords(1:ixyz-1)
    open(eout,file=trim(stem)//'_energy.dat')

    write(eout,'(6F18.7)') 0.0d0,Ekin, Epot,Ekin+Epot
    do i=1,ndt
        if (Nbath > 0) then
            call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, dt, ne)
        end if
        call verletstep(dt, Epot)
        Ekin = kinetic_energy()
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

        write(eout,'(6F18.7)') dt*i*t0fs,Ekin, Epot,Ekin+Epot
    end do
    close(eout)
end program gmd
! vim:et
