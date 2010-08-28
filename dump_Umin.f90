subroutine dump_Umin(MCiter)
    use ljmc
    implicit none
    integer :: MCiter
    integer :: i, j, imin
    
    write (*, "('n =', I10, Umin = F14.7)" MCiter, Umin(1:nstreams)

    imin = minloc(Umin(1:nstreams))

    open(33, file=fname)
    write(33,'(I10)') Natom
    write(33,"('Umin = ', ES14.7)") Umin(imin)
    write(33,"('Ne ', 3ES14.7)") ((r(i, j,imin), i=1,3), j=1,N)
    close(33)
end subroutine

