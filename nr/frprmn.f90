SUBROUTINE frprmn(f, fdf, p, ftol,iter,fret)
    USE nrtype; USE nrutil, ONLY : nrerror
    USE nr, ONLY : linmin, func_s_v, fdf_s_v
    IMPLICIT NONE
    procedure(func_s_v) :: f
    procedure(fdf_s_v) :: fdf
    INTEGER(I4B), INTENT(OUT) :: iter
    REAL(SP), INTENT(IN) :: ftol
    REAL(SP), INTENT(OUT) :: fret
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: p

    INTEGER(I4B), PARAMETER :: ITMAX=200
    REAL(SP), PARAMETER :: EPS=1.0e-10_sp
    INTEGER(I4B) :: its
    REAL(SP) :: dgg,fp,gam,gg
    REAL(SP), DIMENSION(:), allocatable :: g,h,xi

    allocate(g(size(p)), h(size(p)), xi(size(p)))
    call fdf(p, fp, xi)
    g=-xi
    h=g
    xi=h
    do its=1,ITMAX
        iter=its
        call linmin(f, fdf, p,xi,fret)
        if (2.0_sp*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN
!!$        fp=fret
!!$        xi=df(p)
        call fdf(p, fp, xi)
        gg=dot_product(g,g)
        !		dgg=dot_product(xi,xi)
        dgg=dot_product(xi+g,xi)
        if (gg == 0.0) RETURN
        gam=dgg/gg
        g=-xi
        h=g+gam*h
        xi=h
    end do
    call nrerror('frprmn: maximum iterations exceeded')
    deallocate(g,h,xi)
END SUBROUTINE frprmn
