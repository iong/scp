program OHeff
    use vgwfm_mod
    class(vgwfm), pointer :: fm
    double precision :: r0(3,2), Ueff, kT, f(3,2)
    character(256) :: argin
    integer :: ix, nx=100

    double precision, parameter :: A_au = 1.8897, K_au = 3.16681534282952d-06, &
          kOH = 0.49536d0

    if (command_argument_count() >0) then
       call get_command_argument(1, argin)
       read(argin,*) kT
    end if
    if (command_argument_count() >1) then
       call get_command_argument(2, argin)
       read (argin, *) nx
    end if

    kT = kT * K_au

    allocate(fm)
    call fm%init(2,'OH')
    call fm%set_mm(.TRUE.)

    xeq = 0.9824 * A_au !equilibrium of the quantum OH bond
    dx = 3.5d0 /(sqrt(1823.0*kOH)*tanh(0.5d0*sqrt(kOH/1823.0)/kT))**0.5

    dt = 2.0*M_PI*sqrt(1823.0/kOH)/640.0

    open(33,file='OHeff.dat')
    r0 = 0d0
    do ix=1,nx
         r0(1,2) = (2d0*real(ix-1)/real(nx-1) - 1d0)*dx + xeq
         call fm%Ueff(r0, 1d0/kT, Ueff)
         call fm%Feff(f)
         write(33, '(9G14.6)') r0(1,2), Ueff, exp(-Ueff/kT), f
     end do
     close(33)

    write (*,*) dx, xeq
    call fm%cleanup()
end program OHeff
