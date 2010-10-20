program mergecvv
    use utils
    implicit none
    character(len=256) :: fname, cmd, cfgfile, coords
    integer :: Natom, nr0, NMCp, ndt, nr0start, nr0stop, nr0step, i, j, k, Nbath
    logical :: fexists
    real*8 :: M(3,3), Ueff, detgj, Vmin_curvature, Z
    real*8 :: rcmin, tstart, tstop, dtout, dt, bl, bl2,rho, kT
    real*8, dimension(:), allocatable :: y, ekin, epot, etot, Cvv, Cvvavg, Cvvstd, detg, WW, t
    real*8, allocatable :: r0(:,:), G(:,:,:)
    namelist /gmdcfg/Natom,mass,NGAUSS,LJA,LJC,rc,rtol,atol,taumin,kT,rho, &
            rcmin, NMCp,coords,tstart,tstop,dtout,dt,Nbath,Vmin_curvature

    cfgfile='pH2.in'
    taumin = 1.0e-6
    rtol=1e-4
    atol=1.0e-4
    rc=10.0
    LJA=0.0d0
    LJC=0.0d0
    open(30,file=cfgfile)
    read(30,NML=gmdcfg)
    close(30)

    call get_command_argument(1, cmd)
    read(cmd, '(I)') nr0start
    call get_command_argument(2, cmd)
    read(cmd, '(I)') nr0step
    call get_command_argument(3, cmd)
    read(cmd, '(I)') nr0stop

    nr0 = (nr0stop - nr0start) / nr0step + 1
    ndt = (tstop - tstart) / dt  + 1

    allocate ( ekin(0:ndt), epot(0:ndt), etot(0:ndt), Cvv(0:ndt), WW(0:ndt), &
            Cvvavg(0:ndt), Cvvstd(0:ndt), t(0:ndt), detg(nr0))

    Cvvavg = 0.0d0
    Cvvstd = 0.0d0
    Z = 0
    do i=1,nr0
        j = 1
        do
            if (mod(j,2) == 1) then
                write(fname,'("dump/r0_",I10,"_",I5,".dat")') (i-1)*nr0step + nr0start, j+1
                call replace_char(fname, ' ', '0')
                inquire(FILE=fname,EXIST=fexists)
                if (.not. fexists) exit
            end if
            write(fname,'("dump/r0_",I10,"_",I5,".dat")') (i-1)*nr0step + nr0start, j
            call replace_char(fname, ' ', '0')
            open(30,file=fname)
            read(30,*) (t(k), Ekin(k), Epot(k), Etot(k), Cvv(k), k=1,ndt)
            close(30)

            Cvv = Cvv / Natom

            Cvvavg = Cvvavg + Cvv
            Z = Z + 1.0d0
            j = j + 1
        end do
    end do
    Cvvavg = Cvvavg / Z
    
    open(30,file='Cvv.dat')
    write(30,*) (t(i), Cvvavg(i), i=1,ndt)
    close(30)

end program mergecvv
