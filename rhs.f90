module rhs
    interface
        SUBROUTINE RHSS0(NEQ, T, x, xp)
            IMPLICIT NONE
            integer, intent(in) :: NEQ
            REAL*8, intent(in) :: T, x(NEQ)
            REAL*8, intent(out) :: xp(NEQ)
        end subroutine
        SUBROUTINE RHSS1(NEQ, T, x, xp)
            IMPLICIT NONE
            integer, intent(in) :: NEQ
            REAL*8, intent(in) :: T, x(NEQ)
            REAL*8, intent(out) :: xp(NEQ)
        end subroutine
    end interface
end module
