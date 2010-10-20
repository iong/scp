program gmd
    use spine
    use vgw
    use xyz
    use utils
    use propagation
    implicit none
    integer :: NMCp, ndtout, ndt, ne
    real*8 :: v3(3), Ueff0, rmserr, p3(3), Q1nhc, sumf(3), sump(3)
    real*8, allocatable :: f(:,:), WW(:), dr(:,:), dqp(:), pbath0(:)
    real*8 :: Ekin ,Epot, Cvv

    character(LEN=256) :: cfgfile, fname, coords
    integer :: i, j, k, n, ixyz
    namelist /gmdcfg/Natom,mass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,kT,rho, &
            rcmin, NMCp,coords,tstart,tstop,dtout,dt,Nbath,Q1nhc,ne

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
    
    allocate (y(1+21*natom), r0(3,natom), p0(3,natom), vtau0(3,natom), &
                Meff0(3,3,natom), invMeff0(3,3,natom), Qnk0(3,3,natom), &
                Meff(3,3,natom), invMeff(3,3,natom), &
                r(3,natom), p(3,natom), f(3,natom), dr(3,natom))
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
    call vgw1(r0, Ueff0, 1.0/kT, 0.0d0, y, Meff0, invMeff0)
    call unpack_Qnk(y, Qnk0)

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

    do n=1,NMCp
        if (mod(n,2) == 0) then
            p0 = -p0
        else
            call initial_momenta(kT, Meff0, p0)
            !sump =  sum(p0, 2)/Natom
            !do i=1,Natom
            !    p0(:,i) = p0(:,i) - sump
            !end do
        end if
        do i=1,Natom
            v3 = matmul(invMeff0(:,:,i), p0(:,i))
            vtau0(:,i) = matmul(Qnk0(:,:,i), v3)
        end do

        r = r0
        p = p0
        f = 0.0d0

        call total_ekin(Ekin)
        if (Nbath > 0) then
            call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, 0.0d0, ne)
        end if
        call verletstep(r, p, f, Epot, 0.0d0)
        if (Nbath > 0) then
            call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, 0.0d0, ne)
        end if


        ixyz = index(coords, '.xyz', .TRUE.)
        write(fname, "('dump/',A,'_',I5,'.dat')") coords(1:ixyz-1), n
        call replace_char(fname, ' ', '0')
        open(30,file=fname)
        call velocity_autocorrelation(Cvv)
        write(30,'(6F18.7)') 0.0d0, ekin, epot, ekin+epot, Cvv

        do i=1,ndt
            !dqp=qp
            if (Nbath > 0) then
                call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, dt, ne)
            end if

            call verletstep(r, p, f, Epot, dt)
            call total_ekin(Ekin)

            if (Nbath > 0) then
                call nose_hoover_chain(p, Ekin, kT, xi, vxi, Qbath, dt, ne)
                write(31,*) vxi
            end if

            call velocity_autocorrelation(Cvv)
            write(30,'(6F18.7)') dt*i*t0fs,Ekin, Epot,Ekin+Epot, Cvv

            do j=1,Natom
                do k=1,3
                    if (abs(r(k, j)) >bl2) then
                        r(k, j) = r(k, j) - sign(bl, r(k, j))
                    end if
                end do
            end do

            !forall (k=1:3*Natom, abs(qp(k)) >bl2)
            !    qp(k) = qp(k) - sign(bl, qp(k))
            !end forall
        end do
        close(30)
    end do
end program gmd
! vim:et
