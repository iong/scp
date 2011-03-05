program gmdshort
    use spine
    use vgw
    use xyz
    use utils
    use propagation
    use correlations
    implicit none
    integer :: iostat, ip
    integer :: runid=0, np
    real*8 :: Ekin ,Epot, tlen, pcm(3), Ueff

    character(LEN=256) :: cfgfile, coords,argin
    integer :: i, j, k, ixyz
    namelist /gmdcfg/mass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,kT,rho, &
            coords,tstart,tstop,dt,np


    cfgfile='pH2short.in'
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
    if (command_argument_count() >1) then
       call get_command_argument(2, argin)
       read (argin, '(I)') runid 
    end if

    call load_xyz(coords, r0)
    Natom = size(r0, 2)
    ndt = (tstop - tstart) / dt

    allocate (y(1+21*natom), rkold(3, Natom), &
                Qnk(3,3,natom), Meff(3,3,natom), invMeff(3,3,natom), &
                r(3,natom), p(3,natom), p0(3,natom), v(3,natom), f(3, Natom), &
                q0tau(3,natom), rshift(3,natom),r0k(3,natom), &
                v0k(3,natom), r0shift(3,natom), &
                v0tau(3,natom), v0s(3,natom), r0s(3,natom), &
                v0(3,natom), vkubo(3,natom), &
                track(1:track_width,0:ndt), &
                trackaccum(1:track_width,0:ndt))

    call seed_rng()

    bl = (Natom/rho)**(1.0/3.0)
    bl2=bl/2
    mass = mass*0.020614788876D0
    call vgwinit(natom, bl)


    
    ixyz = index(coords, '.xyz', .TRUE.)
    stem = coords(1:ixyz-1)
    open(eout,file='dump/'//trim(stem)//'_energy.dat')

    trackaccum = 0.0d0
    do ip=1,2*np
        r = r0
        call vgw1(r, Ueff, 1.0/kT, 0.0d0, y, Meff, invMeff)
        if (mod(ip,2) == 1) then
            call initial_momenta(kT, Meff, p0)
            pcm = sum(p0, 2)/Natom
            do i=1,Natom
                p0(:,i) = p0(:,i) - pcm
            end do
            p = p0
        else
            p = -p0
        end if
        call verletstep(0.0d0, Epot)
        call init_track()
        call update_track(0)
    
       
        Ekin = kinetic_energy()
        write(eout,'(6F18.7)') 0.0d0,Ekin, Epot,Ekin+Epot
        do i=1,ndt
            call verletstep(dt, Epot)
            call update_track(i)
            call pb_wrap(r, rshift, bl)
            if (ip==1) then
                Ekin = kinetic_energy()
                write(eout,'(6F18.7)') dt*i*t0fs,Ekin, Epot,Ekin+Epot
            end if
        end do
        write (eout, '(//)')
        trackaccum = trackaccum + track
    end do
    trackaccum = trackaccum / (2*Natom*np)
    write (*,*) trackaccum(1,0)
    call dump_track(trackaccum, runid)
    close(eout)
end program gmdshort
! vim:et
