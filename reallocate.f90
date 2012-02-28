module reallocate_mod
    use iso_c_binding

    interface reallocate
        module procedure reallocatei
        module procedure reallocated
    end interface reallocate
contains
    FUNCTION reallocatei(p, n)
        implicit none
        integer, POINTER, DIMENSION(:) :: p, reallocatei
        INTEGER, intent(in) :: n
        INTEGER :: nold, ierr

        ALLOCATE(reallocatei(1:n), STAT=ierr)
        IF(ierr /= 0) STOP "allocate error"
        IF(.NOT. ASSOCIATED(p)) RETURN
        nold = MIN(SIZE(p), n)
        reallocatei(1:nold) = p(1:nold)
        DEALLOCATE(p) 
    END FUNCTION REALLOCATEI

    FUNCTION reallocated(p, n)
        implicit none
        real(c_double), POINTER, DIMENSION(:) :: p, reallocated
        INTEGER, intent(in) :: n
        INTEGER :: nold, ierr

        ALLOCATE(reallocated(1:n), STAT=ierr)
        IF(ierr /= 0) STOP "allocate error"
        IF(.NOT. ASSOCIATED(p)) RETURN
        nold = MIN(SIZE(p), n)
        reallocated(1:nold) = p(1:nold)
        DEALLOCATE(p) 
    END FUNCTION REALLOCATED
end module reallocate_mod
