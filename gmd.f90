program gmd
    use gmd_mod
    use vgw
    use xyz
    use utils
    use propagation
    implicit none
    integer :: NMCp, ndtout, ndt
    real*8 :: v3(3), Ueff0, rmserr, p3(3)
    real*8, allocatable :: r(:,:), v(:,:), a(:,:), WW(:), dr(:,:)

    character(LEN=256) :: cfgfile, fname, coords
    integer :: i, j, k, n
    namelist /gmdcfg/Natom,mass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,kT,rho, &
            rcmin, NMCp,coords,tstart,tstop,dtout,dt

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

    ndt = (tstop - tstart) / dt
    ndtout = dtout/dt
    
    allocate (y(1+21*natom), r0(3,natom), p0(3,natom), vtau0(3,natom),qp(6*natom), &
                ekin(0:ndt), epot(0:ndt), etot(0:ndt), Cvv(0:ndt), WW(0:ndt), &
                Meff0(3,3,natom), invMeff0(3,3,natom), Qnk0(3,3,natom), &
                Meff(3,3,natom), invMeff(3,3,natom), &
                r(3,natom), v(3,natom), a(3,natom), dr(3,natom))

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
    call vgw1(r0, Ueff0, 1.0/kT, 0.0d0, y, Meff, invMeff)
    call unpack_Qnk(y, Qnk0)
    Meff0 = Meff
    invMeff0 = invMeff


    do i=1,Natom
        !Meff0(:,:,i) = 0.1*reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/)*mass, (/3,3/))
        !invMeff0(:,:,i) = 10*reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/)/mass, (/3,3/))
        write (1,*) matmul(Meff0(:,:,i), invMeff0(:,:,i))
    end do
    write (1,*) 'xxx'
    flush(1)

    do n=1,NMCp
        if (mod(n,2) == 0) then
            p0 = -p0
        else
            call initial_momenta(kT, Meff0, p0)
        end if
        do i=1,Natom
            v3 = matmul(invMeff0(:,:,i), p0(:,i))
            vtau0(:,i) = matmul(Qnk0(:,:,i), v3)
            v(:,i) = v3
        end do

        !qp(1:3*Natom) = reshape(r0, (/ 3*Natom /) )
        !qp(1+3*Natom:6*Natom) = reshape(p0, (/ 3*Natom /) )

        r = r0
        call vgw1(r, Ueff0, 1.0/kT, 0.0d0, y, Meff, invMeff)
        do j=1,Natom
            a(:,j)=2.0*kT*matmul(invMeff0(:,:,j), y(2+18*Natom+3*(j-1):1+18*Natom+3*j))
        end do
        !call update_TCF(0)
        WW = 0.0d0
        do i=1,ndt
            !call eulerstep(RHS_Veff, qp, dt, atol, rtol, rmserr)
            v = v + 0.5*dt*a
            dr = v*dt
            r = r + dr
            WW(i) = WW(i-1) + kT*dot_product(reshape(dr, (/3*Natom/)), y(2+18*Natom:1+21*Natom))

            call vgw1(r, Ueff0, 1.0/kT, 0.0d0, y, Meff, invMeff)
            WW(i) = WW(i) + kT*dot_product(reshape(dr, (/3*Natom/)), y(2+18*Natom:1+21*Natom))
            do j=1,Natom
                a(:,j)=2.0*kT*matmul(invMeff0(:,:,j), y(2+18*Natom+3*(j-1):1+18*Natom+3*j))
            end do
            v = v + 0.5*dt*a

            
            epot(i) = -2.0*kT*y(1)
            ekin(i) = 0
            do j=1,Natom
                p3 = matmul(Meff0(:,:,j), v(:,j))
                ekin(i) = ekin(i) + dot_product(p3, v(:,j))
            end do
            ekin(i) = ekin(i)*0.5

            !write (*,*) 'rmserr =', rmserr
            !call update_TCF(i)
            forall (k=1:3*Natom, qp(k) >bl2)
                qp(k) = qp(k) - bl
            end forall
            forall (k=1:3*Natom, qp(k) <-bl2)
                qp(k) = qp(k) + bl
            end forall
        end do

        etot = ekin+epot        

        write(fname, "('dump/traj',I6,'.dat')") n
        call replace_char(fname, ' ', '0')
        open(30,file=fname)
        write(30,'(6F14.7)') (dt*i*t0fs, ekin(i), epot(i), etot(i), WW(i),Cvv(i),i=0,ndt)
        close(30)
        stop
    end do
end program gmd
! vim:et
