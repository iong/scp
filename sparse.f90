module sparse
    use iso_c_binding
    use reallocate_mod
    implicit none
    private

    type,public :: csr
        integer :: nrows, ncols, nnz, nnz_row_max
        integer(c_int), pointer :: ia(:), ja(:), iia(:)
        real(c_double), pointer :: x(:)
    contains
        procedure :: init
        procedure :: reallocate => csr_reallocate
        procedure :: force_symmetry 
        procedure :: mirror_uplo 
        procedure :: from_dense 
        procedure :: dense 
        procedure :: get_3rows_ur 
        procedure :: write_3rows_urwd  
        procedure :: write_3rows_urnd
        procedure :: get_3x3_diag 
        procedure :: set_3x3_diag
        procedure :: gemv
        procedure :: gemm
        procedure :: multiply
        procedure :: multiply_restricted
        procedure :: trace 
        procedure :: find_diagonal 
        procedure :: intersect 
        procedure :: logdet
        procedure :: sort 
        procedure :: cleanup
        procedure :: write => csr_write
    end type csr

    interface
        real(C_DOUBLE) function cholmod_logdet(G, ia, ja, N) BIND(C)
            use, intrinsic :: iso_c_binding
            type(C_PTR), value :: G, ia, ja
            integer(C_INT), value :: N
        end function cholmod_logdet

        subroutine csr_sort_ja(M, ia, ja) bind(c)
            use, intrinsic ::  iso_c_binding
            integer(c_int), value :: M
            type(c_ptr), value :: ia, ja
        end subroutine csr_sort_ja
    end interface

    public :: gemm_restricted

contains

    subroutine init(self, nrows, ncols, nnz, values)
        implicit none
        class(csr), target :: self
        integer, intent(in) :: nrows, ncols, nnz
        logical, intent(in), optional :: values

        logical :: lvalues = .TRUE.

        if (present(values)) then
            lvalues = values
        end if

        self % nrows = nrows
        self % ncols = ncols
        self % nnz = nnz

        allocate(self % ia(nrows + 1), self % iia(nrows), self % ja(nnz) )

        if (lvalues) allocate (self % x (nnz) )

    end subroutine init

    subroutine csr_reallocate(self, n)
        implicit none
        class(csr) :: self
        integer, intent(in) :: n

        self % ja => reallocate(self %ja, n)

        if (associated(self % x)) then
            self % x => reallocate(self % x, n)
        end if
    end subroutine csr_reallocate


    subroutine dense(self, A)
        implicit none
        class(csr) :: self
        double precision, intent(out) :: A(:,:)

        integer :: i, p

        A = 0d0
        do i=1,self%nrows
            do p=self%ia(i) , self%ia(i+1)-1
                A(i, self%ja(p) ) = self%x(p)
            end do
        end do
    end subroutine dense


    subroutine from_dense(self, A)
        implicit none
        class(csr) :: self
        double precision, intent(in) :: A(:,:)

        integer :: i, p

        do i=1,self%nrows

            do p=self%ia(i) , self%ia(i+1)-1

                self%x(p) = A(i, self%ja(p) )

            end do

        end do
    end subroutine from_dense


    subroutine sort(self)
        class(csr), target :: self
        call csr_sort_ja(self % nrows, c_loc(self % ia), c_loc(self % ja) )
    end subroutine sort

    subroutine force_symmetry(self)
        implicit none
        class(csr) :: self

        double precision :: xavg
        integer :: w( self % nrows ), i, j, p, pmirror

        w = 0

        do i=1, self%nrows
            do p=self%ia(i), self%ia(i+1)-1
                j = self%ja(p)
                if (j .le. i) cycle

                pmirror = self%ia(j) + w(j)
                xavg = 0.5d0* (self%x(p) + self%x( pmirror ) )

                self%x(p) = xavg
                self%x( pmirror ) = xavg

                w(j) = w(j) + 1
            end do
        end do
    end subroutine force_symmetry

    subroutine mirror_uplo(self)
        implicit none
        class(csr) :: self

        integer :: w( self % nrows ), i, j, p

        w = 0
        do i=1, self%nrows
            do p=self%iia(i)+1, self%ia(i+1)-1
                j = self%ja(p)

                self % x ( self%ia(j) + w(j) ) = self % x(p)

                w(j) = w(j) + 1
            end do
        end do
    end subroutine mirror_uplo

    !> Get the horizontal 3x band after the diagonal 3x3 block
    ! \param first_row start w this row
    ! \param x matrix data in csr format
    ! \param ia row pointers
    ! \param ja column numbers
    ! \param iia pointers to the diagonal elements
    ! \param band resulting band
    subroutine get_3rows_ur(self, first_row, band)
        implicit none
        class(csr) :: self
        integer, intent(in) :: first_row
        double precision, intent(out) :: band(:,:)

        integer :: rowlen, k, j

        rowlen = self%ia(first_row + 1) - self%ia(first_row)

        k = 0
        do j=self%iia(first_row)+3, self%ia(first_row+1) - 1
            k = k + 1
            band(1, k) = self%x(j)
            band(2, k) = self%x(j+rowlen)
            band(3, k) = self%x(j+2*rowlen)
        end do

        band(:,k+1:) = 0d0

    end subroutine get_3rows_ur



    !> Write the horizontal 3x band including the diagonal 3x3 block
    ! \param first_row start w this row
    ! \param x matrix data in csr format
    ! \param ia row pointers
    ! \param ja column numbers
    ! \param iia pointers to the diagonal elements
    ! \param band resulting band
    subroutine write_3rows_urwd(self, first_row, band)
        implicit none
        class(csr) :: self
        integer, intent(in) :: first_row
        double precision, intent(in) :: band(:,:)

        integer :: ncolcopy

        ncolcopy = self%ia(first_row+1) - self%iia(first_row)

        self%x( self%iia(first_row)       : self%ia(first_row+1) - 1) = band(1,:ncolcopy)
        self%x( self%iia(first_row+1) - 1 : self%ia(first_row+2) - 1) = band(2,:ncolcopy)
        self%x( self%iia(first_row+2) - 2 : self%ia(first_row+3) - 1) = band(3,:ncolcopy)

    end subroutine write_3rows_urwd


    !> Write the horizontal 3x band excluding the diagonal 3x3 block
    ! \param first_row start w this row
    ! \param x matrix data in csr format
    ! \param ia row pointers
    ! \param ja column numbers
    ! \param iia pointers to the diagonal elements
    ! \param band resulting band
    subroutine write_3rows_urnd(self, first_row, band)
        implicit none
        class(csr) :: self
        integer, intent(in) :: first_row
        double precision, intent(in) :: band(:,:)

        integer :: ncolcopy

        ncolcopy = self%ia(first_row+1) - self%iia(first_row) - 3

        self%x( self%iia(first_row) + 3   : self%ia(first_row+1) - 1) = band(1,:ncolcopy)
        self%x( self%iia(first_row+1) + 2 : self%ia(first_row+2) - 1) = band(2,:ncolcopy)
        self%x( self%iia(first_row+2) + 1 : self%ia(first_row+3) - 1) = band(3,:ncolcopy)

    end subroutine write_3rows_urnd


    subroutine find_diagonal(self)
        implicit none
        class(csr) :: self

        integer :: i, p

        do i=1, self%nrows
            do p=self%ia(i), self%ia(i+1) - 1
                if (self%ja(p) == i) then
                    self%iia(i) = p
                endif
            end do
        end do
    end subroutine find_diagonal


    function trace(self)
        implicit none
        class(csr) :: self
        double precision :: trace

        trace = sum(self % x ( self%iia ) )
    end function trace


    !> Get the diagonal 3x3 blocks
    ! \param x matrix data in csr format
    ! \param ia row pointers
    ! \param iia pointers to the diagonal elements
    ! \param diag diagonal 3x3 blocks
    subroutine get_3x3_diag(self, diag)
        implicit none
        class(csr) :: self
        double precision, intent(out) :: diag(:,:)

        integer :: i, N, rowlen, first_col(3)

        do i=1,self%nrows,3
            rowlen = self%ia(i + 1) - self%ia(i)

            first_col = self%iia(i) + (/ 0, 1, 2 /) * rowlen

            diag(:, i)   = self%x(first_col)
            diag(:, i+1) = self%x(first_col + 1)
            diag(:, i+2) = self%x(first_col + 2)
        end do

    end subroutine get_3x3_diag

    !> Set the diagonal 3x3 blocks
    ! \param x matrix data in csr format
    ! \param ia row pointers
    ! \param iia pointers to the diagonal elements
    ! \param diag diagonal 3x3 blocks
    subroutine set_3x3_diag(self, diag)
        implicit none
        class(csr) :: self
        double precision, intent(in) :: diag(:,:)

        integer :: i, N, rowlen, first_col(3)

        do i=1,self%nrows,3
            rowlen = self%ia(i + 1) - self%ia(i)

            first_col = self%iia(i) + (/ 0, 1, 2 /) * rowlen

            self%x(first_col) = diag(:, i) 
            self%x(first_col + 1) = diag(:, i+1)
            self%x(first_col + 2) = diag(:, i+2)
        end do

    end subroutine set_3x3_diag

    subroutine intersect(self, A, idx)
        implicit none
        class(csr) :: self
        type(csr), intent(in) :: A
        integer, intent(out) :: idx(:)

        integer :: i, p, pa

        if (size(idx) < A % nnz) then
            write (*,*) 'Storage is not large enough for intersection!'
            stop
        end if


        do i = 1, self % nrows
            pa = A%ia(i)
            do p = self%ia(i), self%ia(i+1) - 1
                if (self%ja(p) == A%ja (pa) ) then
                    idx(pa) = p
                    pa = pa + 1
                end if
                if (pa >= A%ia(i+1)) exit
            end do
        end do

        if (A%nnz /= pa -1) then
            write (*,*) 'Intersect failed!'
            stop
        end if
    end subroutine intersect

    function logdet(self)
        use iso_c_binding
        implicit none
        class(csr), target :: self
        double precision :: logdet

        logdet = cholmod_logdet(C_LOC(self%x), C_LOC(self%ia), &
                C_LOC(self%ja), self%nrows)
    end function logdet

    subroutine cleanup(self)
        implicit none
        class(csr) :: self
        integer :: istat

        ! possible double deallocation somewhere else
        if (associated(self % x)) then
            deallocate(self%x, STAT=istat)
        end if
        

        if (associated(self%ia)) then
            deallocate (self%ia, self%ja, self%iia)   
        end if

        self % nrows = 0
        self % ncols = 0
        self % nnz = 0
    end subroutine cleanup

    !> y = A*x 
    subroutine gemv (A, x, y)
        implicit none
        class(csr) :: A
        double precision, intent(in) :: x(:)
        double precision, intent(out) :: y(:)

        integer :: i, p

        y = 0d0
        do i = 1, A%nrows
            do p = A%ia(i), A%ia(i+1) - 1
                y(i) = y(i) +  A%x(p) * x( A%ja(p) )
            end do
        end do
    end subroutine gemv

    subroutine gemm (A, x, y)
        implicit none
        class(csr) :: A
        double precision, intent(in) :: x(:,:)
        double precision, intent(out) :: y(:,:)

        integer :: i, p, k

        y = 0d0
        do i = 1, A%nrows
            do k=1,size(x,2)
                do p = A%ia(i), A%ia(i+1) - 1
                    y(i, k) = y(i, k) +  A%x(p) * x( A%ja(p), k )
                end do
            end do
        end do
    end subroutine gemm


    subroutine gemm_restricted(A, B, C, transb, beta)
        double precision, intent(in) :: A(:,:), B(:,:)
        type(csr), intent(inout) :: C
        character(*), intent(in), optional :: transb
        double precision, intent(in), optional :: beta
        integer, parameter :: slw = 12
        double precision, allocatable :: S(:,:)
        double precision :: lbeta=0d0
        character :: ltransb='N'
        integer :: i, j, w, p0, p1, M, N, K

        M = size(A, 1)
        N = size(B, 2)
        K = size(B, 1)

        if (present(transb)) ltransb = transb
        if (present(beta)) lbeta = beta

        if (ltransb == 'T' .or. ltransb == 't') then
            N = size(B, 1)
            K = size(B, 2)
        end if

        allocate(S(slw, N))

        do i=1,M,slw
            w = min(slw, M-i+1)
            call dgemm('N', ltransb, w, N, K, 1d0, A(i,1), size(A, 1), &
                  B, size(B, 1), 0d0, S, slw)
            if (lbeta /= 0d0) then
                do j=1,w
                    do p0 = C%ia(i+j-1), C%ia(i+j) - 1
                        C%x(p0) = lbeta* C%x(p0) + S(j,C%ja(p0))
                    end do
                end do
            else
                do j=1,w
                    p0 = C%ia(i+j-1)
                    p1 = C%ia(i+j) - 1
                    C%x(p0:p1) = S(j,C%ja(p0:p1))
                end do
            end if
        end do

        deallocate(S)
    end subroutine gemm_restricted


    !> x = x + beta * A(i,:), where x is a dense vector and A(:,j) is sparse
    function scatter (A, i, beta, w, x, mark, C, nz)
        implicit none
        class(csr), intent(in) :: A
        integer, intent(in) :: i, mark, nz
        double precision, intent(in) :: beta
        integer, intent(inout) :: w(:)
        double precision, intent(inout) :: x(:)
        integer :: scatter
        type(csr), intent(inout) :: C
        logical :: values

        integer :: Cnz, p, j

        values = ASSOCIATED(C % x) .AND. ASSOCIATED(A % x)

        Cnz = nz
        do p = A % ia(i), A % ia(i+1) - 1
            j = A % ja(p)
            if (w (j) < mark) then
                w (j) = mark
                C%ja(Cnz+1) = j
                Cnz = Cnz + 1
                if (values) then
                    x (j) = beta * A % x (p)
                end if
            elseif (values) then
                x (j) = x(j) + beta * A % x (p)
            end if
        end do
        scatter = Cnz
    end function scatter

    integer function scatter_cnz (A, i, w, mark, nz) result(Cnz)
        implicit none
        class(csr), intent(in) :: A
        integer, intent(in) :: i, mark, nz
        integer, intent(inout) :: w(:)

        integer :: p, j

        Cnz = nz
        do p = A % ia(i), A % ia(i+1) - 1
            j = A % ja(p)
            if (w (j) < mark) then
                w (j) = mark
                Cnz = Cnz + 1
            end if
        end do
    end function scatter_cnz

    function estimate_cnz (A, B) result(Cnz)
        implicit none
        class(csr), intent(in) :: A
        type(csr), intent(in) :: B

        integer, allocatable :: w(:)
        integer :: i, p, Cnz

        allocate (w(A%ncols))
        w = 0

        if (A%ncols /= B%nrows) return

        Cnz = 0
        do i=1,A%nrows
            do p = A%ia(i), A%ia(i+1)-1
                Cnz = scatter_cnz (B, A%ja(p), w, i, Cnz)
            end do
        end do
        deallocate(w)
    end function estimate_cnz


    !>
    !  <pre>
    !  forall i:
    !      forall k:
    !          forall j:
    !              C[i,j] += A[i,k] * B[k,j]
    !  </pre>
    function multiply (A, B) result(C)
        implicit none
        class(csr), intent(in) :: A
        type(csr), intent(in) :: B
        type(csr) :: C

        double precision, allocatable :: x(:)
        integer, allocatable :: w(:)
        logical :: values
        integer :: i, p, Cnz
        
        allocate (w(A%ncols),  x(A%ncols))
        x = 0
        w = 0
        values =  ASSOCIATED(B%x) .AND. ASSOCIATED(A%x)

        if (A%ncols /= B%nrows) return
    
        Cnz = estimate_cnz(A, B)
        call C%init(A%nrows, B%ncols, Cnz, values)

        Cnz = 0
        do i=1,C%nrows
            C%ia(i) = Cnz + 1
            do p = A%ia(i), A%ia(i+1)-1
                Cnz = scatter (B, A%ja(p), A%x(p), w, x, i, C, Cnz)
            end do

            if (values) then
                do p = C % ia(i), Cnz
                    C%x(p) = x ( C%ja (p) )
                end do
            end if
        end do
        C%ia (C%nrows+1) = Cnz + 1
        if (C%nnz /= Cnz) then
            print *, 'multiply error!'
            stop
        end if

        deallocate(w, x)
    end function multiply


 !> x = x + beta * A(i,:), where x is a dense vector and A(:,j) is sparse
    subroutine scatter_restricted (A, i, beta, x)
        implicit none
        class(csr), intent(in) :: A
        integer, intent(in) :: i
        double precision, intent(in) :: beta
        double precision, intent(inout) :: x(:)

        integer :: p, j

        do p = A % ia(i), A % ia(i+1) - 1
            j = A % ja(p)
            x (j) = x(j) + beta * A % x (p)
        end do
    end subroutine scatter_restricted


    !>
    !  <pre>
    !  forall i:
    !      forall k:
    !          forall j:
    !              C[i,j] += A[i,k] * B[k,j]
    !  </pre>
    subroutine multiply_restricted (A, B, C)
        class(csr), intent(in) :: A
        type(csr), intent(in) :: B
        type(csr), intent(inout) :: C

        double precision, allocatable :: x(:)
        integer, allocatable :: w(:)
        integer :: i, p

        

        if (A%ncols /= B%nrows .or. A%nrows /= C%nrows .or. B%ncols /= C%ncols) then
            return
        end if

        C%x = 0d0
!$omp parallel private(w, x)
        allocate (w(A%ncols),  x(A%ncols))
        x = 0
        w = 0
!$omp do private(i, p) schedule(dynamic)
        do i=1,C%nrows
            x = 0d0
            do p = A%ia(i), A%ia(i+1)-1
                call scatter_restricted (B, A%ja(p), A%x(p), x)
            end do

            do p = C % ia(i), C%ia(i+1) - 1
                C%x(p) = x ( C%ja (p) )
            end do
        end do
!$omp end do
        deallocate(w, x)
!$omp end parallel
    end subroutine multiply_restricted


    subroutine csr_write(self, name)
        class(csr) :: self
        character(*) :: name

        open(38, file=trim(name))
        write(38, '(3I10)') self%nrows, self%ncols, self%nnz
        write(38, '(I10)') self%ia, self%ja
        write(38, '(F17.9)'), self%x
        close(38)
    end subroutine csr_write
end module sparse
