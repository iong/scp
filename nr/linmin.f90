SUBROUTINE linmin(f, fdf, p,xi,fret)
    USE nrtype; USE nrutil, ONLY : assert_eq
    USE nr, ONLY : mnbrak,dbrent,func_s_v, fdf_s_v
    IMPLICIT NONE
    procedure(func_s_v) :: f
    procedure(fdf_s_v) :: fdf
    REAL(SP), INTENT(OUT) :: fret
    REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
    REAL(SP), PARAMETER :: TOL=1.0e-4_sp
    REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx

    ax=0.0
    xx=0.05d0
    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
    fret=dbrent(ax,xx,bx,fdf1dim,TOL,xmin)
    xi=xmin*xi
    p=p+xi
contains
    function f1dim(x)
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        REAL(SP) :: f1dim

        f1dim = f(p + xi*x)
    end function f1dim

    subroutine fdf1dim(x, f, df)
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        REAL(SP), INTENT(OUT) :: f, df
        REAL(SP), allocatable :: g(:)

        call fdf(p + xi*x, f, g)
        df = dot_product(g, xi)
        deallocate(g)
    end subroutine fdf1dim
END SUBROUTINE linmin
