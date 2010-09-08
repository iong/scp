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
    character(w), intent(out) :: str0
    integer :: z, n2, i

    n2 = n
    z=iachar('0')
    do i=w,1,-1
        str0(i:i) = achar(z + mod(n2,10))
        n2 = n2/10
    end do
end subroutine

end module utils
