subroutine dump_Umin(MCiter)
    use ljmc
    implicit none
    integer :: MCiter
    integer :: i, j, im(1)
    
    write (*, "('n =',I10,'Umin = ',F14.7)") MCiter, Umin(1:nstreams)

    im = minloc(Umin(1:nstreams))

    open(33, file=outfile)
    write(33,'(I10)') Natom
    write(33,"('Umin = ', ES14.7)") Umin(im(1))
    write(33,"('Ne ', 3ES14.7)") ((r(i, j,im(1)), i=1,3), j=1,Natom)
    close(33)
end subroutine

