module vgwfm_mod
    use integration
    use utils
    use vgw_mod
    implicit none
    private

    type, public, extends(vgw) :: vgwfm
    contains
	procedure :: converge
	procedure :: analyze
	procedure :: logdet
    end type vgwfm
    
    interface
        double precision function dlamch(s)
            character(1) :: s
        end function dlamch
    end interface

contains

    subroutine analyze(self, q0)
        class(vgwfm), target :: self
        double precision, intent(in) :: q0(:)
        double precision :: q

        self % NG =  (9 * self%Natom**2 + 3*self%Natom)/2
        self % NEQ = 3 * self % Natom + self % NG
    end subroutine analyze

    subroutine converge(self, tol)
        IMPLICIT NONE
        class(vgwfm) :: self
        double precision, intent(in) :: tol
        double precision, allocatable :: Omega(:,:), Hnew(:,:), G(:,:)

        integer :: N3

        N3 = 3 * self % Natom

        allocate(Omega(N3,N3), Hnew(N3,N3), G(N3, N3))
        
        call self % prop % converge(RHS, self%y, tol)

        deallocate(Omega, Hnew, G)
    contains

        SUBROUTINE rhs(NEQ, T, Y, YP)
            !    use omp_lib
            IMPLICIT NONE
            integer, intent(in) :: NEQ
            double precision, intent(in) :: T
            double precision, intent(in), target :: Y(:)
            double precision, intent(out), target :: YP(:)

            integer :: i, j

            if (y(N3+1) == 0d0) then
                yp = 0d0
                j = 1
                Hnew = 0d0
                do i=1,N3
                    Hnew(i,i) = 1d0
                end do
                !all regtransrot(y(1:N3), Hnew, 0d0)
                call fm_set_g(Hnew, yp)
                return
            end if

            call fm_get_g(y, Omega)
            
            !call printev(Omega, 'Omega')
            call dsyrk('L', 'N', N3, N3, 2d0 * self % kT, Omega, N3, 0d0, G, N3)
            call GaussianAverage(self, y(1:N3), G, self%U, self % UX, Hnew)   

            call dsymv('L', 'N', -1d0, Omega, N3, self % UX, 1, 0d0, yp, 1)

            call dsymm('L', 'L', N3, N3, 1d0, Omega, N3, Hnew, N3, 0d0, G, N3)
            Hnew = 0d0
            do i=1,N3
                Hnew(i,i) = 1d0
            end do
            call dsymm('R', 'L', N3, N3, -1d0, Omega, N3, G, N3, 1d0, Hnew, N3)

            call regtransrot(y(1:N3), Hnew, 0d0)
            call fm_set_g(Hnew, yp)

            self % qconv  = sqrt(sum(yp(1:N3)**2)/N3)/(2d0*self%kT)
            self % gconv = sqrt(sum(yp(N3+1:)**2)/(NEQ-N3))

            !print *, T, self % qconv, self %gconv 
        END SUBROUTINE RHS
    end subroutine CONVERGE


    SUBROUTINE GaussianAverage(self, q, G, U, UX, UXY)
        !    use omp_lib
        IMPLICIT NONE
        type(vgwfm) :: self
        double precision, intent(in) :: q(:), G(:,:)
        double precision, intent(out), target :: U, UX(:), UXY(:,:)
        INTEGER :: J,I1,I2,IG, N3
        REAL*8 AG(3,3), &
                DETA,DETAG,QZQ,U12, &
                G12(3,3),A(3,3), &
                Zq(3), Z(3,3),Q12(3)
        real*8 :: UXY0(3,3), UX0(3)

        N3 = size(q)

        U = 0; UX = 0; UXY = 0;
        do I1 = 1, N3 - 3,3
            do I2 = I1+3, N3, 3
                Q12 = q(I1:I1+2) - q(I2:I2+2)
                Q12 = min_image(Q12, self%bl)
                G12=  g_cross_block(G, I1, I2)

                call detminvm(G12, DETA, A)
                DETA = 1.0d0/DETA

                UX0 = 0d0; UXY0 = 0d0
                DO IG=1, self%NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
                    AG = A
                    do J=1,3
                        AG(J,J) = self%LJA(IG) + AG(J,J)
                    end do

                    call detminvm(AG, DETAG, Z)
                    Z = - self%LJA(IG)**2 * Z

                    do J=1,3
                        Z(J,J) = Z(J,J) + self%LJA(IG)
                    end do

                    Zq = matmul(Z, Q12) ! R = -2.0*Zq
                    qZq = dot_product(Q12, Zq)  

                    if (DETA/DETAG <0) then
                        print *, 'ltzero'
                    endif
                    U12 = SQRT(DETA/DETAG) * EXP(-qZq) * self%LJC(IG)
                    U = U + U12

                    UX0 = UX0 - 2d0*U12*Zq
                    do J=1,3
                        UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                    end do
                end do ! IG

                ! Avoid load and store as much as possbile. Store now,
                ! process as a stream later. Much faster.
                UX (I1 : I1 + 2) = UX (I1 : I1 + 2) + UX0
                UX (I2 : I2 + 2) = UX (I2 : I2 + 2) - UX0

                UXY0 = 0.5d0 * (UXY0 + transpose(UXY0))

                UXY(I1 : I1+2, I1 : I1+2) = UXY(I1 : I1+2, I1 : I1+2) + UXY0
                UXY(I2 : I2+2, I2 : I2+2) = UXY(I2 : I2+2, I2 : I2+2) + UXY0

                UXY (I1 : I1+2, I2 : I2+2) = -UXY0
                UXY (I2 : I2+2, I1 : I1+2) = -UXY0
            end do ! I2
        end do ! I1
    END SUBROUTINE GaussianAverage

    subroutine fm_get_g(y, G)
        implicit none
        double precision, intent(in) :: y(:)
        double precision, intent(out) :: G(:,:)

        integer :: ypos, i, N

        N = size(G, 1)

        ypos = N + 1
        do i = 1, N
            G( i:, i) = y (ypos : ypos + N - i)
            G( i, i : N) = y (ypos : ypos + N - i)
            ypos = ypos + N - i + 1
        end do
    end subroutine fm_get_g

    subroutine fm_set_g(G, y)
        implicit none
        double precision, intent(in) :: G(:,:)
        double precision, intent(out) :: y(:)

        integer :: ypos, i, N

        N = size(G, 1)

        ypos = N + 1
        do i = 1,N
            y (ypos : ypos + N - i) = G(i : N, i)
            ypos = ypos + N - i + 1
        end do
    end subroutine fm_set_g

    function logdet(self)
        class(vgwfm), target :: self
        double precision :: logdet

        double precision, allocatable :: G(:,:)

        integer :: info, j

        allocate(G(3 * self%Natom, 3 * self%Natom))
        call fm_get_g(self % y, G)
        call regtransrot(self%y(1:3*self%Natom), G, 1d0)

        call dpotrf('U', 3 * self%Natom, G, 3 * self%Natom, info)

        logdet=0.0
        DO j=1, 3 * self%Natom
            logdet = logdet + LOG(ABS( G(j,j) ))
        ENDDO
        logdet = 2d0* logdet
        deallocate(G)
    end function logdet


    subroutine regtransrot(q, G, l)
        implicit none
        double precision :: q(:), G(:,:), l

        double precision, allocatable :: U(:,:), GU(:,:), UGU(:,:)
        integer :: i, j, N

        N = size(G, 1)

        allocate(U(N,6), GU(N,6), UGU(6,6))

        U = transrot_subspace(q)

        call dgemm('N', 'N', N, 6, N, 1d0, G, N, U, N, 0d0, GU, N)
        call dgemm('T', 'N', 6, 6, N, -1d0, U, N, GU, N, 0d0, UGU, 6)

        do i=1,6
            UGU(i,i) = UGU(i,i) + l
        end do

        !print '(6F12.7)', UGU

        call dgemm('N', 'N', N, 6, 6, 1d0, U, N, UGU, 6, 0d0, GU, N)
        call dgemm('N', 'T', N, N, 6, 1d0, GU, N, U, N, 1d0, G, N)
        deallocate(U, GU, UGU)
    end subroutine regtransrot

    pure function g_cross_block(G, I, J) result(x)
        double precision, intent(in) :: G(:,:)
        integer, intent(in) :: I, J
        double precision ::x(3,3)
        x = G(I : I+2, I : I+2) + G(J : J+2, J : J+2)
        x(1,2) = x(2,1)
        x(1,3) = x(3,1)
        x(2,3) = x(3,2)
        if (I>J) then
            x = x - G(I : I+2, J : J+2) - transpose(G(I : I+2, J : J+2)) 
        else
            x = x - G(J : J+2, I : I+2) - transpose(G(J : J+2, I : I+2))
        end if
    end function g_cross_block

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
end module vgwfm_mod
