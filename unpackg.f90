subroutine unpackg(N_atom, y, G)
    integer, intent(in) :: N_atom    
    real*8, intent(in) :: y(6*N_atom)
    real*8, intent(out) :: G(3,3,N_atom)
    integer :: i, j, k, cnt

    cnt = 1
    do i=1,N_atom
        do j=1,3
            do k=j,3
                G(k, j, i) = y(cnt)
                G(j, k, i) = y(cnt)
		cnt = cnt + 1
	    enddo
	enddo
    enddo
end subroutine
