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
    real*8 :: x(2)
    do
        call random_number(x)
        if (x(1) /= 0.0d0) exit
    enddo
    !   write (*,*) x

    y = sqrt( -2.0 * log(x(1))) * cos(2*M_PI*x(2))

    y = y*sigma + x0
end function gaussran

subroutine pol2cart(pol, cart)
    real*8, intent(in) :: pol(3)
    real*8, intent(out) :: cart(3)
    real*8 :: rxy

    rxy = pol(1) * sin(pol(2))
    cart(1) = rxy * cos(pol(3))
    cart(2) = rxy * sin(pol(3))
    cart(3) = pol(1) * cos(pol(2))
end subroutine pol2cart


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
end subroutine int2strz

subroutine linspace(xmin, xmax, N, xout)
    real*8, intent(in) :: xmin, xmax
    integer, intent(in) :: N
    real*8, intent(out) :: xout(N)
    integer :: i
    real*8:: dx

    dx = (xmax - xmin) / (N-1)
    xout = xmin + dx * (/(i,i=0,N-1)/)
end subroutine linspace

subroutine replace_char(str, a, b)
    character(len=*), intent(inout) :: str
    character, intent(in) :: a, b
    integer :: i

    do i=1,len_trim(str)
        if (str(i:i)==a) then
            str(i:i) = b
        end if
    end do
end subroutine replace_char

subroutine fliplr(v)
    real*8, intent(inout) :: v(:)
    real*8 :: x
    integer :: i, N

    N = size(v)
    do i=1,N/2
        x = v(N-i+1)
        v(N-i+1) = v(i)
        v(i) = x
    end do
end subroutine fliplr

subroutine detminvm(A, DETA, INVA)
    implicit none
    real*8, intent(in) :: A(3,3)
    real*8, intent(out) :: DETA, INVA(3,3)

    INVA(1,1) = A(2,2)*A(3,3)-A(2,3)*A(3,2)
    INVA(2,1) = -A(1,2)*A(3,3)+A(1,3)*A(3,2)
    INVA(3,1) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
    INVA(1,2) = INVA(2,1)
    INVA(2,2) = A(1,1)*A(3,3)-A(1,3)*A(3,1)
    INVA(3,2) = -A(1,1)*A(2,3)+A(1,3)*A(2,1)
    INVA(1,3) = INVA(3,1)
    INVA(2,3) = INVA(3,2)
    INVA(3,3) = A(1,1)*A(2,2)-A(1,2)*A(2,1)

    DETA = INVA(1,1)*A(1,1)+INVA(2,1)*A(2,1)+INVA(3,1)*A(3,1)
    INVA = INVA / DETA
end subroutine detminvm

function outer_product(l, r) result(m)
    real*8, intent(in) :: l(:), r(:)
    real*8 :: m(size(l), size(r))
    integer :: i

    forall (i=1:size(r))
        m(:,i) = l*r(i)
    end forall
end function outer_product

pure function outer_product3(l, r) result(m)
    real*8, intent(in) :: l(3), r(3)
    real*8 :: m(3,3)
    integer :: i

    forall (i=1:3)
        m(:,i) = l*r(i)
    end forall
end function outer_product3

pure subroutine pb_wrap(r, rshift, bl)
    real*8, intent(in) :: bl
    real*8, intent(inout) :: r(:,:), rshift(:,:)
    integer :: j, k
    real*8 :: bl2

    bl2 = 0.5d0 * bl
    do j=1,size(r, 2)
        do k=1,3
            if (abs(r(k, j)) > bl2) then
                rshift(k, j) = rshift(k, j) + sign(bl, r(k, j))
                r(k, j) = r(k, j) - sign(bl, r(k, j))
            end if
        end do
    end do
end subroutine pb_wrap
end module utils
