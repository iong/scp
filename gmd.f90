program gmd
    use spine
    use vgw
    use xyz
    use utils
    use propagation
    implicit none
    integer :: iostat
    integer :: ndtout, ndt, ne, nequil
    real*8 :: v3(3), Ueff0, rmserr, p3(3), Q1nhc, sumf(3), sump(3)
    real*8, allocatable :: f(:,:)
    real*8 :: Ekin ,Epot, Cvv, tequil, dxsq,tlen, pcm(3)

    character(LEN=256) :: cfgfile, coords
    integer :: i, j, k, n, ixyz
    namelist /gmdcfg/Natom,mass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,kT,rho, &
            rcmin,coords,tstart,tstop,dtout,dt,Nbath,Q1nhc,ne,tequil,tlen


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
    ndtout = dtout/dt
    nequil = tequil / dt
    window_width = tlen / dt

    naccumulated = 0
    ncvvout = 0
    
    allocate (y(1+21*natom), r0(3,natom), &
                Qnk(3,3,natom), Meff(3,3,natom), invMeff(3,3,natom), &
                r(3,natom), p(3,natom), f(3,natom), v(3,natom),&
                r0equil(3,natom),rshift(3,natom),&
                v0tau(3,natom,window_width), v0s(3,natom,window_width), &
                p0(3,natom,window_width), &
                Cvvcur(3,window_width), Cvvold(3,window_width))
    call load_xyz(r0, coords)
    call seed_rng()

    bl = (Natom/rho)**(1.0/3.0)
    bl2=bl/2

    if (rcmin**3 * rho >= 1.0) then
        write (*,"(A,F7.3)") 'rcmin is too large. Ensure that rcmin <', &
            rho**(-1.0/3.0)
        stop
    endif

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
    call verletstep(r0, p, f, Epot, 0.0d0) ! compute Meff
    call initial_momenta(kT, Meff, p)
    pcm = sum(p, 2)/Natom
    do i=1,Natom
        p(:,i) = p(:,i) - pcm
    end do

    r = r0
    f = 0.0d0
    v0tau = 0.0d0
    Cvvold = 0.0d0
    Cvvcur = 0.0d0

    call total_ekin(Ekin)
    if (Nbath > 0) then
        call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, 0.0d0, ne)
    end if
    call verletstep(r, p, f, Epot, 0.0d0)
    if (Nbath > 0) then
        call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, 0.0d0, ne)
    end if


    ixyz = index(coords, '.xyz', .TRUE.)
    write (*,*) coords(1:ixyz-1), ixyz
    stem = 'dump/'//coords(1:ixyz-1)
    write (*,*) stem
    open(eout,file=trim(stem)//'_energy.dat')

    write(eout,'(6F18.7)') 0.0d0,Ekin, Epot,Ekin+Epot
    do i=1,ndt
        if (Nbath > 0) then
            call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, dt, ne)
        end if
        call verletstep(r, p, f, Epot, dt)
        call total_ekin(Ekin)
        if (Nbath > 0) then
            call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, dt, ne)
        end if

        if (i>=nequil) then
            call velocity_autocorrelation()
        end if
        write(eout,'(6F18.7)') dt*i*t0fs,Ekin, Epot,Ekin+Epot


        do j=1,Natom
            do k=1,3
                if (abs(r(k, j)) >bl2) then
                    r(k, j) = r(k, j) - sign(bl, r(k, j))
                    rshift(k, j) = rshift(k, j) + sign(bl, r(k, j))
                end if
            end do
        end do
    end do
    close(eout)
end program gmd
! vim:et
