module mainvars
    use vgwfm_mod
    implicit none
    logical :: debug = .FALSE.

    real*8, parameter :: t0 = 7.6382d-12, t0fs = 7638.2d0
    integer, parameter :: eout = 30, cvvout = 31, mapout=32, track_width = 14
    integer :: Natom, ndt
    real*8 :: rcmin, tstart, tstop, dt, bl, bl2,rho, kT

    real*8, dimension(:), allocatable :: mass
    real*8, dimension(:,:), allocatable :: r0, r, p0, p, v, rshift, f, &
        trackaccum, rkold, v0tau, v0s, v0, vkubo, track, q0tau, r0shift, &
        r0k, rprev, vprev, vsprev, r0s, rsprev,v0k
    real*8, dimension(:,:), allocatable :: Qnk, Meff, invMeff, sqrtMeff, &
        sqrtInvMeff

    real*8 :: lastepot
    character(256) :: stem


    type(vgwfm), allocatable, target :: fm
end module mainvars
