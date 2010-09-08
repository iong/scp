module xyz
contains
subroutine dump_xyz(r, fname, molname)
    real*8, intent(in) :: r(:,:)
    character(256), intent(in) :: fname, molname
    integer :: i, j, N

    N = size(r, 2)
    open(33,file=fname)
    write(33,'(I10)') N
    write(33,*) trim(molname)
    write(33,"(('Ne ', 3(F12.5, ' ')))") ((r(i, j), i=1,3), j=1,N)
    close(33)
end subroutine

subroutine load_xyz(r, fname)
    real*8, intent(out) :: r(:,:)
    character(256), intent(in) :: fname
    character(256):: molname
    integer :: i, j, N

    open(33,file=fname)
    read(33,*) N
    if (N /= size(r, 2)) then
        write (*,*) "N /= size(r, 2)", N, size(r, 2)
        stop
    end if
    read(33,"(A)") molname
    read(33, *) (molname, r(:, j), j=1,N)
    close(33)
end subroutine
end module xyz
