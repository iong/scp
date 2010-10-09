program gmd
    use spine
    use vgw
    use xyz
    use utils
    use propagation
    implicit none
    integer :: NMCp, ndtout, ndt
    real*8 :: v3(3), Ueff0, rmserr, p3(3), Vmin_curvature
    real*8, allocatable :: r(:,:), p(:,:), f(:,:), WW(:), dr(:,:), dqp(:), pnh0(:)

    character(LEN=256) :: cfgfile, fname, coords
    integer :: i, j, k, n, ixyz
    namelist /gmdcfg/Natom,mass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,kT,rho, &
            rcmin, NMCp,coords,tstart,tstop,dtout,dt,Nbath,Vmin_curvature

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
    
    allocate (y(1+21*natom), r0(3,natom), p0(3,natom), vtau0(3,natom),qp(6*natom+Nbath), dqp(6*natom+Nbath),&
                ekin(0:ndt), epot(0:ndt), etot(0:ndt), Cvv(0:ndt), WW(0:ndt), &
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

    if (Nbath>0) then
        allocate(Qbath(Nbath), pnh0(Nbath))
        Qbath = kT*mass/Vmin_curvature
        Qbath(1) = 3*Natom*Qbath(1)
        do i=1,Nbath
            pnh0(i) = gaussran(sqrt(kT*Qbath(i)), 0.0d0)
        end do
    end if

    do n=1,NMCp
        if (mod(n,2) == 0) then
            p0 = -p0
        else
            call initial_momenta(kT, Meff0, p0)
        end if
        do i=1,Natom
            v3 = matmul(invMeff0(:,:,i), p0(:,i))
            vtau0(:,i) = matmul(Qnk0(:,:,i), v3)
        end do

        qp(1:3*Natom) = reshape(r0, (/ 3*Natom /) )
        qp(1+3*Natom:6*Natom) = reshape(p0, (/ 3*Natom /) )
        if (Nbath>0) then
            qp(1+6*Natom:6*Natom+Nbath) = pnh0
        end if

        call update_TCF(0)
        WW = 0.0d0

        ixyz = index(coords, '.xyz', .TRUE.)
        write(fname, "('dump/',A,'_',I5,'.dat')") coords(1:ixyz-1), n
        call replace_char(fname, ' ', '0')
        open(30,file=fname)

        etot = 0.d0
        epot = 0.0d0
        ekin = 0.0d0
        Cvv = 0.0d0
        do i=1,ndt
            dqp=qp
            call eulerstep(RHS_Veff, qp, dt, atol, rtol, rmserr)
            dqp = qp - dqp
            WW(i) = WW(i-1) + 2.0*kT*dot_product(dqp(1:3*Natom), y(2+18*Natom:1+21*Natom))

            write (*,*) 'rmserr =', rmserr
            call update_TCF(i)
            write(30,'(6F18.7)') dt*i*t0fs, ekin(i), epot(i), etot(i), WW(i),Cvv(i)

            forall (k=1:3*Natom, qp(k) >bl2)
                qp(k) = qp(k) - bl
            end forall
            forall (k=1:3*Natom, qp(k) <-bl2)
                qp(k) = qp(k) + bl
            end forall
        end do

        
        !write(30,'(6F14.7)') (dt*i*t0fs, ekin(i), epot(i), etot(i), WW(i),Cvv(i),i=0,ndt)
        close(30)
    end do
end program gmd
! vim:et
