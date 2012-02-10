MODULE nr
    abstract INTERFACE
        FUNCTION func_s_s(p)
            USE nrtype
            IMPLICIT NONE
            REAL(SP), INTENT(IN) :: p
            REAL(SP) :: func_s_s
        END FUNCTION func_s_s

        FUNCTION func_s_v(p)
            USE nrtype
            IMPLICIT NONE
            REAL(SP), DIMENSION(:), INTENT(IN) :: p
            REAL(SP) :: func_s_v
        END FUNCTION func_s_v

!!$        FUNCTION func_v_s(p)
!!$            USE nrtype
!!$            IMPLICIT NONE
!!$            REAL(SP), INTENT(IN) :: p
!!$            REAL(SP), DIMENSION(size(p)) :: dfunc
!!$        END FUNCTION func_v_s

        FUNCTION func_v_v(p)
            USE nrtype
            IMPLICIT NONE
            REAL(SP), DIMENSION(:), INTENT(IN) :: p
            REAL(SP), DIMENSION(size(p)) :: func_v_v
        END FUNCTION func_v_v

        subroutine fdf_s_s(x, f, df)
            USE nrtype
            IMPLICIT NONE
            REAL(SP), INTENT(IN) :: x
            REAL(SP), INTENT(OUT) :: f, df
        end subroutine fdf_s_s

        subroutine fdf_s_v(x, f, df)
            USE nrtype
            IMPLICIT NONE
            REAL(SP), INTENT(IN) :: x(:)
            REAL(SP), INTENT(OUT) :: f
            REAL(SP), INTENT(OUT), allocatable :: df(:)
        end subroutine fdf_s_v

        FUNCTION brent(ax,bx,cx,func,tol,xmin)
            USE nrtype
            import func_s_s
            REAL(SP), INTENT(IN) :: ax,bx,cx,tol
            REAL(SP), INTENT(OUT) :: xmin
            REAL(SP) :: brent
            procedure(func_s_s) :: func
        END FUNCTION brent

        FUNCTION dbrent(ax,bx,cx,fdf,tol,xmin)
            USE nrtype
            IMPORT fdf_s_s
            REAL(SP), INTENT(IN) :: ax,bx,cx,tol
            REAL(SP), INTENT(OUT) :: xmin
            REAL(SP) :: dbrent
            procedure(fdf_s_s) :: fdf
        END FUNCTION dbrent


        SUBROUTINE frprmn(f, fdf, p,ftol,iter,fret)
            USE nrtype
            IMPORT FUNC_S_V
            IMPORT fdf_s_v
            procedure(func_s_v) :: f
            procedure(fdf_s_v) :: fdf
            INTEGER(I4B), INTENT(OUT) :: iter
            REAL(SP), INTENT(IN) :: ftol
            REAL(SP), INTENT(OUT) :: fret
            REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
        END SUBROUTINE frprmn

        SUBROUTINE linmin(f, fdf, p,xi,fret)
            USE nrtype
            import func_s_v
            import fdf_s_v
            procedure(func_s_v) :: f
            procedure(fdf_s_v) :: fdf
            REAL(SP), INTENT(OUT) :: fret
            REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
        END SUBROUTINE linmin

        SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
            USE nrtype
            IMPORT FUNC_S_S
            REAL(SP), INTENT(INOUT) :: ax,bx
            REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
            procedure(func_s_s) :: func
        END SUBROUTINE mnbrak
    END INTERFACE
END MODULE nr
