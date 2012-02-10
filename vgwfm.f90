module vgwfm_mod
    use integration
    use utils
    use vgw_mod
    implicit none
    private

    type, public, extends(vgw) :: vgwfm
	real*8 :: rfullmatsq
	real*8, allocatable :: UXY(:,:)
    contains
	procedure :: cleanup
	procedure :: propagate
	procedure :: init_prop
	procedure :: logdet
	procedure :: hessian
    end type vgwfm
    
    interface
        double precision function dlamch(s)
            character(1) :: s
        end function dlamch
    end interface

contains

    subroutine cleanup(self)
        implicit none
        class(vgwfm) :: self

        deallocate(self%UXY, self%y)
        
        call self%vgw%cleanup()
    end subroutine cleanup


    subroutine init_prop(self, q0)
        class(vgwfm), target :: self
        double precision, intent(in) :: q0(:,:)
        integer :: N3

        self % NG =  9 * self%Natom**2
        self % NEQ = 3 * self % Natom! + self % NG

        N3 = 3 * self % Natom
        if (.NOT. allocated(self % y)) then
            allocate(self % y ( self % NEQ ), self%UXY(N3, N3))
        end if

        allocate (lsode :: self%prop)

        call self % prop % init(self%NEQ, 0d0, self%dt0)
        call self % prop % set_dtmin(0d0*self%dt_min)
        call self % prop % set_dtmax(1d-1)

        self % prop % rtol = 1d-5
        self % prop % atol (1 : 3 * self % Natom) = 1d-3!self % q_atol
        !self % prop % atol (3 * self % Natom + 1 : ) = 1d-5!self % g_atol

        self % y(1 : N3) = reshape(Q0, (/ N3 /) )

        !G(1:N3,1:N3) => self%y(N3+1:)

        call self % hessian( self % y(1 :N3), self%UXY)
        !call add_external_field( self % y(1 :N3), self%UXY)
        self%UXY = 0.5d0 * (self%UXY + transpose(self%UXY))
    end subroutine init_prop

    subroutine propagate(self, tstop)
        IMPLICIT NONE
        class(vgwfm) :: self
        double precision, intent(in) :: tstop
        
        call self % prop % advance(RHS, self%y, tstop)
    contains

        SUBROUTINE rhs(NEQ, T, Y, YP)
            !    use omp_lib
            IMPLICIT NONE
            integer, intent(in) :: NEQ!, IPAR(:)
            double precision, intent(in) :: T
            double precision, intent(in), target :: Y(:)!, RPAR(:)
            double precision, intent(out), target :: YP(:)

            double precision :: dUXYnorm
            double precision, allocatable :: UXYnew(:,:)

            allocate(UXYnew,source=self%UXY)
 
            do
                call avg_ref_system(self, y, self%U, yp, UXYnew)
                dUXYnorm = sqrt(sum((UXYnew - self % UXY)**2)/size(UXYnew))
                !print *, 'dUXYnorm = ', dUXYnorm
                if (dUXYnorm < 1d-2) then
                    exit
                end if
                self % UXY = UXYnew
            end do
            
            yp=-yp

            self % qconv = sqrt(sum(yp**2)/size(yp))

            print *,  T, self % qconv!, self % gconv

            deallocate(UXYnew)
        END SUBROUTINE RHS
    end subroutine PROPAGATE

    SUBROUTINE avg_ref_system(self, q, U, UX, UXY)
        !    use omp_lib
        IMPLICIT NONE
        type(vgwfm) :: self
        double precision, intent(in) :: q(:)
        double precision, intent(out), target :: U, UX(:), UXY(:,:)
        INTEGER :: J,I1,I2,IG,  N3, info
        REAL*8 AG(3,3), &
                DETA,DETAG,QZQ,U12, &
                G12(3,3),A(3,3), &
                Zq(3), Z(3,3),Q12(3)
        real*8 :: UXY0(3,3), UX0(3)

        double precision, allocatable :: G(:,:)

        N3 = 3*self%Natom

        allocate(G(N3,N3))
        G = pinv(self%UXY)
!!$        G = self%UXY
!!$        call dpotrf('L', N3, G, N3, info)
!!$        call dpotri('L', N3, G, N3, info)
!!$
!!$        do i1=1,N3
!!$            G(i1,i1:) = G(i1:,i1)
!!$        end do

        U = 0; UX = 0; UXY = 0;
        do I1 = 1, 3*(self%Natom - 1),3
            do I2 = I1+3, 3*self%Natom, 3
                Q12 = q(I1:I1+2) - q(I2:I2+2)
                Q12 = min_image(Q12, self%bl)
                G12=  G (I1 : I1+2, I1 : I1+2) + G (I2 : I2+2, I2 : I2+2) &
                        - G (I2 : I2+2, I1 : I1+2) - G (I1 : I1+2, I2 : I2+2)
                G12 = 2d0 * self%kT * G12

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
                    self % U = self % U + U12

                    UX0 = UX0 - 2d0*U12*Zq
                    do J=1,3
                        UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                    end do
                end do ! IG

                ! Avoid load and store as much as possbile. Store now,
                ! process as a stream later. Much faster.
                UX (I1 : I1 + 2) = UX (I1 : I1 + 2) + UX0
                UX (I2 : I2 + 2) = UX (I2 : I2 + 2) - UX0

                 UXY(I1 : I1 + 2, I1 : I1 + 2) = &
                         UXY(I1 : I1 + 2, I1 : I1 + 2) + UXY0

                UXY (I2 : I2 + 2, I2 : I2 + 2) = &
                        UXY (I2 : I2 + 2, I2 : I2 + 2) + UXY0

                UXY (I1 : I1 + 2, I2 : I2 + 2) = -UXY0
                UXY (I2 : I2 + 2, I1 : I1 + 2) = -transpose(UXY0)
            end do ! I2
        end do ! I1

        !call add_external_field( self % y(1 :N3), UXY)
        deallocate(G)
    END SUBROUTINE AVG_REF_SYSTEM

    function pinv(A)
        double precision :: A(:,:)
        double precision :: pinv(size(A,2),size(A,1))
        double precision, allocatable :: AA(:,:), work(:), W(:), U(:,:)
        integer,allocatable :: iwork(:),isuppz(:)

        integer :: N, i, info, NEV

        N = size(A, 1)

        allocate(AA,source=A)
        allocate(work(N * 26), W(N), U(N,N), iwork(N * 10), isuppz(2*N))
        call dsyevr('V', 'A', 'L', N, AA, N, 0d0, 0d0, 1, 6, &
                DLAMCH('Safe minimum'), NEV, W, U, N, isuppz, &
                work, size(work), iwork, size(iwork), info)

        AA(:,1:6) = 0d0
        forall (i=7:N) AA(:,i) = U(:,i) / W(i)
        call dgemm('N', 'T', N, N, N, 1d0, U, N, AA, N, 0d0, pinv, N)

        deallocate(AA,W,U,work,iwork,isuppz)
    end function pinv

    function pinvs(A)
        double precision :: A(:,:)
        double precision :: pinvs(size(A,2),size(A,1))
        double precision, allocatable :: AA(:,:), work(:), W(:), U(:,:)
        integer,allocatable :: iwork(:),isuppz(:)

        integer :: N, i, info, NEV

        N = size(A, 1)

        allocate(AA,source=A)
        allocate(work(N * 26), W(N), U(N,N), iwork(N * 10), isuppz(2*N))
        call dsyevr('V', 'A', 'L', N, AA, N, 0d0, 0d0, N-5, N, &
                DLAMCH('Safe minimum'), NEV, W, U, N, isuppz, &
                work, size(work), iwork, size(iwork), info)

        AA(:,N-5:N) = 0d0
        forall (i=1:N-6) AA(:,i) = U(:,i) / W(i)
        call dgemm('N', 'T', N, N, N, 1d0, U, N, AA, N, 0d0, pinvs, N)

        deallocate(AA,W,U,work,iwork,isuppz)
    end function pinvs

    subroutine printev(A)
        double precision :: A(:,:)
        double precision, allocatable :: AA(:,:), work(:), W(:)
        integer,allocatable :: iwork(:),isuppz(:)

        integer :: N, info, NEV
        interface
            double precision function dlamch(s)
                character(1) :: s
            end function dlamch
        end interface

        N = size(A, 1)

        allocate(AA,source=A)
        allocate(work(N * 26), W(N), iwork(N * 10), isuppz(2*N))
        call dsyevr('N', 'A', 'L', N, AA, N, 0d0, 0d0, 1, 6, &
                DLAMCH('Safe minimum'), NEV, W, 0, N, isuppz, &
                work, size(work), iwork, size(iwork), info)
        print '(15F14.6)',W(1:9),W(N-5:N)

        deallocate(AA,W,work,iwork,isuppz)
    end subroutine printev

    subroutine slash_g(q, G)
        double precision :: q(:), G(:,:)

        double precision, allocatable :: U(:,:), GU(:,:)
        integer :: i, j, N

        N = size(G, 1)

        print *, 'slash!!'

        allocate(U(N,6), GU(N,6))
        U = 0d0
        U(1::3,1) = 1d0/sqrt(N/3d0)
        U(2::3,2) = 1d0/sqrt(N/3d0)
        U(3::3,3) = 1d0/sqrt(N/3d0)
        U(2::3,4) = -q(3::3)
        U(3::3,4) =  q(2::3)
        U(1::3,5) =  q(3::3)
        U(3::3,5) = -q(1::3)
        U(1::3,6) = -q(2::3)
        U(2::3,6) =  q(1::3)

        do i=4,6
            do j=1,i-1
                U(:,i) = U(:,i) - U(:,j) * sum ( U(:,i) * U(:,j) )
            end do
            U(:,i) = U(:,i) / sqrt(sum(U(:,i)**2))
        end do

        call dsymm('L', 'L', N, 6, 1d0, G, N, U, N, 0d0, GU, N)
        call dgemm('N', 'T', N, N, 6, -1d0, U, N, GU, N, 1d0, G, N)
        deallocate(U, GU)
    end subroutine slash_g

    subroutine slash_g2(G, which)
        double precision :: G(:,:)
        character(1), intent(in), optional :: which

        double precision, allocatable :: U(:,:), GU(:,:)
        integer :: il, iu, N

        double precision, allocatable :: AA(:,:), work(:), W(:)
        integer,allocatable :: iwork(:),isuppz(:)

        integer :: info, NEV

        N = size(G, 1)

        print *, 'slash2!!'

        allocate(AA,source=G)
        allocate(U(N,6), GU(N,6), work(N * 26), W(N), iwork(N * 10), &
                isuppz(2*N))

        il = 1; iu=6
        if (present(which)) then
            if (which == 'L') then
                il = N-5; iu = N
            end if
        end if
        call dsyevr('V', 'I', 'L', N, AA, N, 0d0, 0d0, il, iu, &
                DLAMCH('Safe minimum'), NEV, W, U, N, isuppz, &
                work, size(work), iwork, size(iwork), info)


        call dsymm('L', 'L', N, 6, 1d0, G, N, U, N, 0d0, GU, N)
        call dgemm('N', 'T', N, N, 6, -1d0, U, N, GU, N, 1d0, G, N)
        deallocate(U, GU, AA,W,work,iwork,isuppz)
    end subroutine slash_g2

                

    subroutine hessian(self, q, UXY)
        !    use omp_lib
        IMPLICIT NONE
        class(vgwfm) :: self
        double precision, intent(in) :: q(:)
        double precision, intent(out) :: UXY(:,:)

        double precision  :: U12, Q12(3), UXY0(3,3)
        integer :: J,I1,I2,IG

        do I1=1,size(q),3
            DO I2=I1+3,size(q),3
                Q12 = min_image(q(I1 : I1+2) - q(I2 : I2+2), self%bl)

                UXY0 = 0d0
                DO IG=1,self%NGAUSS
                    U12 = EXP(-self%LJA(IG) * sum(Q12**2) ) * self%LJC(IG)

                    do J=1,3
                        UXY0(:,J) = UXY0(:,J) + 4d0*U12*self%LJA(IG)**2 * Q12(j) * Q12
                    end do

                    do J=1,3
                        UXY0(J,J) = UXY0(J,J) - 2d0*U12 * self%LJA(IG)
                    end do
                end do


                UXY(I1 : I1+2, I1 : I1+2) = UXY(I1 : I1+2, I1 : I1+2) + UXY0
                UXY(I2 : I2+2, I2 : I2+2) = UXY(I2 : I2+2, I2 : I2+2) + UXY0

                UXY(I1 : I1+2, I2 : I2+2) = -UXY0
                UXY(I2 : I2+2, I1 : I1+2) = -UXY0
            end do ! I2
        end do ! I1
    END SUBROUTINE hessian

    subroutine add_external_field(q, UXY)
        double precision, intent(in) :: q(:)
        double precision, intent(inout) :: UXY(:,:)

        double precision, allocatable :: U(:,:)
        integer :: i, j, N

!!$        do i=1,3
!!$            UXY(i,i) = UXY(i,i) + 0.1d0
!!$        end do

        N = size(q)

        allocate(U(N,6))

        U = 0d0
        U(1::3,1) = 1d0/sqrt(N/3d0)
        U(2::3,2) = 1d0/sqrt(N/3d0)
        U(3::3,3) = 1d0/sqrt(N/3d0)
        U(2::3,4) = -q(3::3)
        U(3::3,4) =  q(2::3)
        U(1::3,5) =  q(3::3)
        U(3::3,5) = -q(1::3)
        U(1::3,6) = -q(2::3)
        U(2::3,6) =  q(1::3)

        do i=4,6
            do j=1,i-1
                U(:,i) = U(:,i) - U(:,j) * sum ( U(:,i) * U(:,j) )
            end do
            U(:,i) = U(:,i) / sqrt(sum(U(:,i)**2))
        end do

        call dgemm('N', 'T', N, N, 6, .1d0, U, N, U, N, 1d0, UXY, N)
        deallocate(U)
    end subroutine add_external_field
        

    subroutine fm_get_g(N, y, G)
        implicit none
        integer, intent(in) :: N
        double precision, intent(in) :: y(:)
        double precision, intent(out) :: G(:,:)

        integer :: ypos, i

        ypos = 3 * N + 1
        do i = 1, 3 * N
            G( i : 3 * N, i) = y (ypos : ypos + 3*N - i)
            G( i, i : 3 * N) = y (ypos : ypos + 3*N - i)
            ypos = ypos + 3 * N - i + 1
        end do
    end subroutine fm_get_g

    subroutine fm_set_g(N, G, y)
        implicit none
        integer, intent(in) :: N
        double precision, intent(in) :: G(:,:)
        double precision, intent(out) :: y(:)

        integer :: ypos, i

        ypos = 3 * N + 1
        do i = 1,3 * N
            y (ypos : ypos + 3 * N - i) = G(i : 3 * N, i)
            ypos = ypos + 3 * N - i + 1
        end do
    end subroutine fm_set_g

    
    function logdet(self)
        class(vgwfm), target :: self
        double precision :: logdet

        double precision, allocatable :: G(:,:)

        integer :: info, j

        allocate(G,source=self%UXY)
        call dpotrf('U', 3 * self%Natom, G, 3 * self%Natom, info)

        logdet=0.0
        DO j=1, 3 * self%Natom
            logdet = logdet + LOG(ABS( G(j,j) ))
        ENDDO
        logdet = 1d0/ (2d0* logdet)
    end function logdet

end module vgwfm_mod
