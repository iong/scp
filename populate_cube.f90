subroutine populate_cube(bl, rcmin, r)
    real*8, intent(in) :: bl, rcmin
    real*8, intent(out) :: r(:,:)
    real*8 :: dsq, rcminsq
    integer :: i, j
    logical :: too_close

    rcminsq = rcmin*rcmin
    do i=1,size(r,2)
        do
            too_close = .FALSE.
            call random_number(r(:,i))
            r(:,i) = (r(:,i) - 0.5) * bl
            do j=1,i-1
                dsq = sum((r(:,j)-r(:,i))**2)
                if (dsq < rcminsq) then
                    too_close = .TRUE.
                    exit
                endif
            enddo
            if (.not. too_close) exit
        enddo
    enddo
end subroutine
