program gmd
    use xyz
    use utils
   ! use propagation
    implicit none
    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer :: NMCp, ndtout, ndt, Nbath, Natom, Ngauss
    real*8 :: ekin, epot, Cvv, sump(3), mass, bl, bl2, rho, rtol, atol, rc, LJA(20), LJC(20), taumin, tstart, tstop, dt, dtout, rcmin, kT,Vmin_curvature
    real*8, allocatable :: r0(:,:), p0(:,:), r(:,:), p(:,:), f(:,:)
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
    
    allocate (r0(3,natom), p0(3,natom), r(3,natom), p(3,natom), f(3,natom))
    call load_xyz(r0, coords)
    call seed_rng()

    bl = (Natom/rho)**(1.0/3.0)
    bl2=bl/2
    mass = mass*0.020614788876D0
    ndt = (tstop - tstart) / dt

    do n=1,NMCp
        if (mod(n,2) == 0) then
            p0 = -p0
        else
            do i=1,Natom
                do j=1,3
                    p0(j,i) = gaussran(sqrt(kT*mass), 0.0d0)
                end do
            end do
            sump =  sum(p0, 2)/Natom
            do i=1,Natom
                p0(:,i) = p0(:,i) - sump
            end do
        end if

        r = r0
        p=p0

        ixyz = index(coords, '.xyz', .TRUE.)
        write(fname, "('dump/lj_',A,'_',I5,'.dat')") coords(1:ixyz-1), n
        call replace_char(fname, ' ', '0')
        open(30,file=fname)

        call lj_f(r, f, epot)
        ekin = 0.5*sum(sum(p0**2, 2))/mass
        Cvv = sum(sum(p*p0, 2))/mass**2
        write(30,'(5F18.7)') 0.0d0, ekin, epot, ekin+epot,Cvv

        do i=1,ndt
            
            p = p + 0.5*dt*f
            r = r + (dt/mass)*p

            call lj_f(r, f, epot)
            p = p + 0.5*dt*f

            ekin = 0.5*sum(sum(p**2, 2))/mass
            Cvv = sum(sum(p*p0, 2))/mass**2

            write(30,'(5F18.7)') dt*i*t0fs, ekin, epot, ekin+epot, Cvv
            write(31,*) reshape(r, (/ 3*Natom /))

            do j=1,Natom
                do k=1,3
                    if (abs(r(k, j)) >bl2) then
                        r(k, j) = r(k, j) - sign(bl, r(k, j))
                    end if
                end do
            end do
        end do
        close(30)
    end do
    stop
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

end program gmd
! vim:et
