subroutine dump_Tmin(MCiter)
    use ljmc
    use mc
    implicit none
    integer :: MCiter
    integer :: i, j
    
    write (*, "('n =',I10,', Umin = ',F14.7)") MCiter, Umin!(1:nstreams)
    write (*,*)


    open(33, file=outfile)
    write(33,'(I10)') Natom
    write(33,"('Umin = ', ES14.7)") Umin
    write(33,"(('Ne ', 3(F12.5, ' ')))") ((rmin(i, j), i=1,3), j=1,Natom)
    close(33)
end subroutine

