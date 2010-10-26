program mergecvv
    use utils
    implicit none
    integer :: iostat
    character(len=256) :: fname, cmd, cfgfile, coords
    integer :: nr0, ndt, nr0start, nr0stop, nr0step, i, j, k,nequil
    logical :: fexists
    real*8 :: tstart, tstop, dt,tequil,Z
    real*8, dimension(:), allocatable :: ekin, epot, etot, Cvv, Cvvavg, Cvvstd, t
    real*8, allocatable :: r0(:,:), dxsq(:), dxsqavg(:)
    namelist /gmdcfg/tstart,tstop,dt,tequil

    open(30,file='pH2.in')
    read(30,NML=gmdcfg,IOSTAT=iostat)
    close(30)

    call get_command_argument(1, cmd)
    read(cmd, '(I)') nr0start
    call get_command_argument(2, cmd)
    read(cmd, '(I)') nr0step
    call get_command_argument(3, cmd)
    read(cmd, '(I)') nr0stop

    nr0 = (nr0stop - nr0start) / nr0step + 1
    ndt = (tstop - tstart) / dt
    nequil = tequil / dt

    allocate ( ekin(0:ndt), epot(0:ndt), etot(0:ndt), Cvv(0:ndt), &
            Cvvavg(0:ndt), Cvvstd(0:ndt), t(0:ndt), dxsq(0:ndt), dxsqavg(0:ndt))

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
            read(30,*) (t(k), Ekin(k), Epot(k), Etot(k), Cvv(k), dxsq(k), k=1,ndt-nequil+1)
            close(30)

            Cvvavg = Cvvavg + Cvv
            dxsqavg = dxsqavg + dxsq
            Z = Z + 1.0d0
            j = j + 1
        end do
    end do
    Cvvavg = Cvvavg / Z
    dxsqavg = dxsqavg / Z
    
    open(30,file='Cvv.dat')
    write(30,'(3F14.7)') (t(i), Cvvavg(i), dxsqavg(i), i=1,ndt-nequil+1)
    close(30)

end program mergecvv
