SUBROUTINE linmin(f, df, p,xi,fret)
    USE nrtype; USE nrutil, ONLY : assert_eq
    USE nr, ONLY : mnbrak,dbrent,func_s_v, func_v_v
    IMPLICIT NONE
    procedure(func_s_v) :: f
    procedure(func_v_v) :: df
    REAL(SP), INTENT(OUT) :: fret
    REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
    REAL(SP), PARAMETER :: TOL=1.0e-4_sp
    REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx

    ax=0.0
    xx=1.0
    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
    fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
    xi=xmin*xi
    p=p+xi
contains
    function f1dim(x)
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        REAL(SP) :: f1dim

        f1dim = f(p + xi*x)
    end function f1dim

    function df1dim(x)
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        REAL(SP) :: df1dim

        df1dim = dot_product(df(p + xi*x), xi)
    end function df1dim
END SUBROUTINE linmin
