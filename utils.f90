module utils
    real*8 :: M_PI = 3.14159265358979323846264338327950288d0
contains
subroutine seed_rng()
        integer :: seed(128), seed_size, c, crate, cmax

        call random_seed(SIZE=seed_size)
!        write(*,*) 'Seed size =', seed_size
        call random_seed(GET=seed(1:seed_size))
        call system_clock(c, crate, cmax)
        seed(1) = mod(c,crate)
        call random_seed(PUT=seed(1:seed_size))
end subroutine

real*8 function gaussran(sigma, x0) result(y)
    implicit none
    real*8, intent(in) :: sigma, x0
    real*8 :: x(2), M_PI = 3.14159265358979323846264338327950288
    call random_number(x)
    do while (x(1) == 0)
        call random_number(x)
    enddo
    y = sqrt( -2.0 * log(x(1))) * cos(2*M_PI*x(2))

    y = y*sigma + x0
end function

subroutine int2strz(n, w, str0)
    implicit none
    integer, intent(in) :: n, w
    character(len=*), intent(out) :: str0
    integer :: z, n2, i

    if (len(str0) < w) then
        write (*,*) 'Number too large for string'
        stop
    end if

    n2 = n
    z=iachar('0')
    do i=w,1,-1
        str0(i:i) = achar(z + mod(n2,10))
        n2 = n2/10
    end do
end subroutine

subroutine linspace(xmin, xmax, N, xout)
    real*8, intent(in) :: xmin, xmax
    integer, intent(in) :: N
    real*8, intent(out) :: xout(N)
    integer :: i
    real*8:: dx

    dx = (xmax - xmin) / (N-1)
    xout = xmin + dx * (/(i,i=0,N-1)/)
end subroutine

subroutine replace_char(str, a, b)
    character(len=*), intent(inout) :: str
    character, intent(in) :: a, b
    integer :: i

    do i=1,len_trim(str)
        if (str(i:i)==a) then
            str(i:i) = b
        end if
    end do
end subroutine

end module utils
