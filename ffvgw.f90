module ffvgw_mod
    use utils
    use sparse
    use vgw_mod
    implicit none
    private

    type, public, extends(vgw) :: ffvgw
        double precision :: Vcutoff, Gcutoff
        type(csr) :: Omega, UXY, OmegaU
    contains
        procedure :: cleanup
        procedure :: rhs
        procedure :: set_range
        procedure :: analyze
        procedure :: logdet
    end type ffvgw

    type csr_list
        type(csr), pointer :: p
    end type csr_list

contains

    subroutine set_range(self, Gcutoff, Vcutoff)
        implicit none
        class(ffvgw) :: self
        double precision, intent(in) :: Vcutoff, Gcutoff

        self % Vcutoff = Vcutoff
        self % Gcutoff = Gcutoff
    end subroutine set_range

    subroutine cleanup(self)
        class(ffvgw) :: self
        call self % Omega   % cleanup()
        call self % UXY % cleanup()
        call self % OmegaU  % cleanup()

        call self % vgw % cleanup()
    end subroutine cleanup


    function count_pairs(Q, srad, bl) result(npairs)
        implicit none
        real*8, intent(in) :: Q(:), srad(:), bl
        integer ::npairs(size(srad))
        integer :: N,I,J, NN
        real*8 :: rsq,sradsq(size(srad)),qij(3)

        N = size(Q)
        sradsq=srad**2

        npairs = 0
        do I=1,N-3,3
            do J=I+3,N,3
                qij=Q(I:I+2)-Q(J:J+2)
                rsq = sum(min_image(qij, BL)**2)
                where (sradsq >= rsq)
                    npairs = npairs + 1
                end where
            enddo
        enddo
    end function count_pairs


    !> Init the sparse pattern
    !  \param ia row pointer
    !  \param ja column number
    subroutine init_sparse_patterns(self, Q0, srad, ML)
        class(ffvgw)::self
        double precision, intent(in) :: Q0(:), srad(:)
        type(csr_list) :: ML(:)
        integer :: I, J, k, ib, jb, Nrad
        integer, dimension(size(srad)) :: p, p0, nnz_row
        double precision :: qij(3), rsq, sradsq(size(srad))

        sradsq = srad**2
        p = 1
        Nrad = size(srad)
        do k=1,Nrad
            ML(k) % p % nnz_row_max = 0
        end do
        do ib=1,self%Natom
            i = 3*ib - 2

            p0 = p

            do jb=1,self%Natom
                j = 3*jb - 2

                qij=Q0(i:i+2)-Q0(j:j+2)
                rsq = sum(min_image(qij, self%BL)**2)

                do k=1,Nrad
                    if (rsq <= sradsq(k) ) then
                        if (ib == jb) then
                            ML(k)%p %iia(i) = p(k)
                        end if

                        ML(k)%p %ja(p(k) : p(k)+2) = (/ j, j+1, j+2/)
                        p(k) = p(k) + 3
                    end if
                end do
            enddo

            nnz_row = p - p0

            do k=1, Nrad
                ML(k)%p % nnz_row_max = max(ML(k)%p % nnz_row_max, nnz_row(k) )
                !write (*,*) ML%p % nnz_row_max, nnz_row
                ML(k)%p % ia(i : i+2) = (/ p0(k), p(k), p(k) + nnz_row(k) /)

                ML(k)%p % ja(p(k) : p(k) + nnz_row(k) - 1) = ML(k)%p % ja(p0(k) : p(k) - 1) 
                ML(k)%p % ja(p(k) + nnz_row(k) : p(k) + 2*nnz_row(k) - 1) = ML(k)%p % ja(p0(k) : p(k) - 1) 

                ML(k)%p % iia(i+1 : i+2) = ML(k)%p % iia(i) + (nnz_row(k) + 1) * (/1, 2/)
            end do

            p = p + 2*nnz_row
        enddo

        do k=1,Nrad
            ML(k)%p % ia(3*self%Natom + 1) = p(k)

            if (ML(k)%p % nnz /= p(k) - 1) then
                write (*,*) 'nnz error!', k, ML(k)%p%nnz, p(k)-1
                stop
            end if
        end do
    end subroutine init_sparse_patterns

    function Utot0(self, Q0) result(U)
        implicit none
        class(ffvgw) :: self
        double precision, intent(in) :: Q0(:,:)
        double precision :: U
        INTEGER  I,J,N
        real*8 :: rsq, QIJ(3)

        N = size(Q0, 2)
        
        U=0d0
        DO I=1,N-1
            DO J=I+1,N
                qij = Q0(:,I) - Q0(:,J)
                rsq = sum(min_image(qij, self % BL)**2)

                if (rsq > self%Vcutoff**2) cycle

                U = U + sum( self % LJC (1 : self % NGAUSS) * &
                        EXP(-self % LJA ( 1 : self % NGAUSS ) * rsq) )
            ENDDO
        ENDDO
    end function Utot0


    SUBROUTINE analyze(self, Q0)
        IMPLICIT NONE
        class(ffvgw), target :: self
        double precision, intent(in) :: Q0(:)

        integer :: nnz(2), N3

        type(csr_list) :: matrix_list(2)


        call self % vgw % analyze(q0)

        N3 = 3 * self % Natom

        nnz = count_pairs(q0, (/ self%Vcutoff, self%Gcutoff /), self%bl )
        nnz = 18*nnz + 3*N3

        call self % UXY % cleanup()
        call self % Omega % cleanup()
        call self % OmegaU % cleanup()

        call self % UXY%init(N3, N3, nnz(1))
        call self % Omega%init(N3, N3, nnz(2))

        matrix_list(1)%p => self%UXY
        matrix_list(2)%p => self%Omega
        call init_sparse_patterns(self, Q0, (/ self%Vcutoff, self%Gcutoff /), matrix_list )

        self % OmegaU = self % Omega %multiply(self % UXY)
        call self % OmegaU % sort()

        if (debug >= info) then
            print '("G  sparsity = ",F7.3,"%")', 100d0 * self % Omega % nnz &
                  / ( self % Omega%nrows * self % Omega%ncols)
            print '("U  sparsity = ",F7.3,"%")', 100d0 * self % UXY%nnz &
                  / ( self % UXY%nrows * self % UXY%ncols)
            print '("GU sparsity = ",F7.3,"%")', 100d0 * self % OmegaU%nnz &
                  / ( self % OmegaU % nrows * self % OmegaU % ncols)
        end if

        deallocate(self % Omega%x)

        self % NEQ = N3 + self % Omega%nnz      
    end subroutine analyze


    subroutine test_ia_ja(ia, ja)
        integer, intent(in) :: ia(:), ja(:)
        integer,allocatable :: P(:,:)
        integer :: i, ptr, N

        allocate(P(size(ia)-1, size(ia)-1))

        N = size(ia) - 1

        P = 0
        do i=1,N
            do ptr=ia(i), ia(i+1)-1
                P(i, ja(ptr)) = 1
            end do
        end do

        if (sum(P-transpose(P)) /= 0) then
            write (*,*) 'The pattern is not symmetric!'
            stop
        end if

        deallocate(P)
    end  subroutine test_ia_ja

    function rhs(self, Y, T)
        IMPLICIT NONE
        class(ffvgw) :: self
        double precision, intent(in) :: T
        double precision, intent(in), target :: Y(:)
        double precision, target :: rhs(size(y))

        type(csr) :: OmegaUOmega, G

        integer :: i, N3

        N3 = 3 * self%Natom

        OmegaUOmega = self % Omega
        OmegaUOmega%x => rhs(N3+1 : N3 + self%Omega%nnz)

        if (y(N3+1)==0d0) then
            rhs = 0d0
            OmegaUOmega%x(OmegaUOmega %iia) = 1d0
            !call regtransrot(y(1:N3), OmegaUOmega, 0d0)
            return
        end if

        self%Omega%x => y(N3+1 : N3 + self%Omega%nnz)

        G = self % Omega
        allocate(G%x(G%nnz))
!$omp parallel
        call self%Omega%multiply_restricted(self%Omega, G)
!$omp do schedule(static)
        do i=1,G%nnz
            G%x(i) = G%x(i) * 2d0 * self % kT
        end do
!$omp end do
        call GaussianAverage(self, y, G, self%U, self%UX, self%UXY)


        !call self% Omega % gemv(self%UX, rhs, 2d0 * self % kT)
        call self% Omega % gemv(self%UX, rhs, -1d0)
        call self%Omega % multiply_restricted(self%UXY, self%OmegaU)
        call self%OmegaU % multiply_restricted(self%Omega, OmegaUOmega)
!$omp end parallel
        OmegaUOmega%x = -OmegaUOmega%x
        OmegaUOmega%x (OmegaUOmega%iia) = OmegaUOmega%x (OmegaUOmega%iia) + 1d0

        call OmegaUOmega % force_symmetry()
        !call regtransrot(y(1:N3), OmegaUOmega, 0d0)
        nullify(self%Omega%x)

        self % qconv = sqrt(sum(rhs(1:N3)**2)/ self % Natom)
        self % gconv = sqrt(sum(rhs(N3:)**2)/ OmegaUOmega % nnz)

        deallocate(G%x)
    end function rhs

    subroutine work_range(UXY, r0, r1)
        use omp_lib
        implicit none
        type(csr) :: UXY
        integer, intent(out) :: r0, r1
        integer :: nthreads, n, i, wrow, wsum
        double precision :: wavg
        integer, allocatable :: row(:)
        
        nthreads = omp_get_num_threads()
        allocate(row(nthreads+1))

        wavg = real(sum(UXY%iia(::3) - UXY%ia(1:UXY%nrows : 3))) / nthreads
        wavg = wavg / 3

        n = 1
        row(n) = 1
        wsum = 0
        do i=1,UXY%nrows,3
            wrow = (UXY%ia(i+1) - UXY%iia(i))/3
            if (abs(wsum - wavg*n) < abs(wsum + wrow - wavg*n) ) then
                n = n+1
                row(n) = i
            end if
            wsum = wsum + wrow
        end do
        row(nthreads+1) = UXY%nrows-2 ! skip the last particle!

        i = omp_get_thread_num()
        r0 = row(i+1)
        r1 = row(i+2)
        deallocate(row)
    end subroutine


    SUBROUTINE GaussianAverage(self, q, G, U, UX, UXY)
        IMPLICIT NONE
        class(ffvgw) :: self
        double precision, intent(in) :: q(:)
        type(csr), intent(in) :: G
        double precision, intent(out) :: U, UX(3 *self%Natom)
        type(csr), intent(out) :: UXY

        INTEGER :: J,I1,I2,IG, p, p1, p2, N3
        double precision :: AG(3,3), DETA,DETAG,QZQ,U12, Gb12(3,3),A(3,3), &
              Zq(3), Z(3,3),Q12(3), UXY0(3,3), UX0(3)

        integer :: Gbandpos

        double precision, allocatable :: Gband(:,:), Gbdiag(:,:)
        integer, allocatable :: Gband_ja(:)

        N3 = 3 * self%Natom
        allocate(Gbdiag(3,N3), Gband(3,G%nnz_row_max), Gband_ja(G%nnz_row_max))

        call G%get_3x3_diag(Gbdiag)

!$omp single
        U = 0; UX = 0;
!$omp end single
!$omp do schedule(dynamic) reduction(+:U,UX)
        do I1=1,N3-3, 3
            Gbandpos = 1
            Gband_ja = 0
            if (G%iia(I1) + 3 < G%ia(I1+1)) then
                call G%get_3rows_ur(I1, Gband)
                Gband_ja(1 : G%ia(I1+1) - G%iia(I1) - 3) = &
                      G%ja(G%iia(I1) + 3 : G%ia(I1 + 1) - 1)
            end if

            do p=UXY%iia(I1) + 3, UXY%ia(I1+1)-1, 3 
                I2 = UXY%ja(p)

                Q12 = q(I1 : I1 + 2) - q(I2 : I2 + 2)
                Q12 = min_image(Q12, self % bl)

                Gb12=Gbdiag(:,I1 : I1+2) + Gbdiag(:,I2 : I2+2)
                if ( Gbandpos < G%nnz_row_max .and. Gband_ja(Gbandpos) == I2) then

                    Gb12 = Gb12 - Gband(:,Gbandpos : Gbandpos + 2) &
                          - transpose(Gband(:,Gbandpos : Gbandpos + 2))
                    Gbandpos = Gbandpos + 3
                end if

                call detminvm(Gb12, DETA, A)
                DETA = 1.0d0/DETA

                UX0 = 0d0; UXY0 = 0d0
                DO IG=1,self%NGAUSS
                    AG = A
                    do J=1,3
                        AG(J,J)=self%LJA(IG)+AG(J,J)
                    end do

                    call detminvm(AG, DETAG, Z)
                    Z = - self%LJA(IG)**2 * Z

                    do J=1,3
                        Z(J,J) = Z(J,J) + self%LJA(IG)
                    end do

                    Zq = matmul(Z, Q12) ! R = -2.0*Zq
                    qZq = dot_product(Q12, Zq) 

                    if (DETA*DETAG <= 0.0d0 ) then
                        print *, 'ltzero'
                        stop
                    end if
                    U12 = SQRT(DETA/DETAG)*EXP(-qZq)*self%LJC(IG)
                    U = U + U12

                    UX0 = UX0 - 2d0*U12*Zq
                    do J=1,3
                        UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                    end do
                end do ! IG

                UX(I1 : I1 + 2) = UX(I1 : I1 + 2) + UX0
                UX(I2 : I2 + 2) = UX(I2 : I2 + 2) - UX0

                UXY0 = -UXY0
                call UXY%put_3x3_block(I1, p, UXY0)
            end do ! I2

        end do ! I1
!$omp end do
        call UXY%mirror_uplo()
        deallocate (Gbdiag, Gband_ja, Gband )

!$omp do schedule(static)
        do i1=1,N3
            p = UXY%iia(i1) - mod(i1-1,3)
            p1 = UXY%ia(i1)
            p2 = UXY%ia(i1+1) - 1
            UXY%x(p : p+2) = 0d0

            UX0(1) = -sum(UXY%x(p1 : p2 : 3))
            UX0(2) = -sum(UXY%x(p1+1 : p2 : 3))
            UX0(3) = -sum(UXY%x(p1+2 : p2 : 3))

            UXY%x(p : p+2) = UX0
        end do
!$omp end do
    end subroutine GaussianAverage


    function logdet(self)
        implicit none
        class(ffvgw), target :: self
        double precision :: logdet
        
        integer :: N3

        N3 = 3 * self%Natom

        self % Omega % x => self % y(N3 + 1 : N3 + self % Omega % nnz)
  
        logdet =  self % Omega %logdet() &
                - logdet_transrot(self%y(1:N3), self % Omega)

        nullify( self % Omega % x)
    end function logdet


    subroutine regtransrot(q, G, l)
        implicit none
        double precision :: q(:), l
        type(csr), intent(inout) :: G

        double precision, allocatable :: U(:,:), GU(:,:), UGU(:,:)
        integer :: i, j, N

        N = G%nrows

        allocate(U(N,6), GU(N,6), UGU(6,6))

        U = transrot_subspace(q)
        call G%gemm(U, GU)
        call dgemm('T', 'N', 6, 6, N, -1d0, U, N, GU, N, 0d0, UGU, 6)

        do i=1,6
            UGU(i,i) = UGU(i,i) + l
        end do

        !print '(6F12.7)', UGU

        call dgemm('N', 'N', N, 6, 6, 1d0, U, N, UGU, 6, 0d0, GU, N)

        call gemm_restricted(GU, U, G, 'T', 1d0)
        deallocate(U, GU, UGU)
    end subroutine regtransrot

    function logdet_transrot(q, G)
        implicit none
        double precision :: q(:)
        type(csr), intent(in) :: G
        double precision :: logdet_transrot

        double precision, allocatable :: U(:,:), GU(:,:), UGU(:,:)
        integer :: j, N, info

        N = G%nrows

        allocate(U(N,6), GU(N,6), UGU(6,6))

        U = transrot_subspace(q)
        call G%gemm(U, GU)
        call dgemm('T', 'N', 6, 6, N, 1d0, U, N, GU, N, 0d0, UGU, 6)

        call dpotrf('U',6, UGU, 6, info)

        logdet_transrot=0.0
        DO j=1, 6
            logdet_transrot = logdet_transrot + LOG(ABS( UGU(j,j) ))
        ENDDO
        logdet_transrot = 2d0 * logdet_transrot
        deallocate(U, GU, UGU)
    end function logdet_transrot

     subroutine printev(A, name)
        double precision :: A(:,:)
        character(*) :: name
        double precision, allocatable :: AA(:,:), work(:), W(:), U(:,:)
        integer,allocatable :: iwork(:),isuppz(:)

        integer :: N, info, NEV
        interface
            double precision function dlamch(s)
                character(1) :: s
            end function dlamch
        end interface

        N = size(A, 1)

        allocate(AA,source=A)
        allocate(U(N,N))
        allocate(work(N * 26), W(N), iwork(N * 10), isuppz(2*N))
        call dsyevr('V', 'A', 'L', N, AA, N, 0d0, 0d0, 1, 6, &
                DLAMCH('Safe minimum'), NEV, W, U, N, isuppz, &
                work, size(work), iwork, size(iwork), info)
        print '(A, 15F14.6)',name, W(1:9),W(N-5:N)

        deallocate(AA,W,work,iwork,isuppz, U)
    end subroutine printev
end module ffvgw_mod
