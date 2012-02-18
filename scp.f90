module scp_mod
    use integration
    use nr, only: frprmn
    use utils, only : M_PI, detminvm

    implicit none
    private

    integer :: N3, Natom

    double precision :: kT
    double precision, pointer :: qcur(:), G(:,:), Glin(:), Hnew(:,:), GH(:,:), &
            LJC(:), LJA(:)

    interface
        double precision function dlamch(s)
            character(1) :: s
        end function dlamch
    end interface


    public:: scp
contains

    subroutine scp(r0, LJC_, LJA_, kT_, F, iter)
        implicit none
        double precision, intent(inout) :: r0(:)
        double precision, intent(in) :: LJC_(:), LJA_(:), kT_
        double precision, intent(out) :: F

        integer, intent(out) :: iter

        class(integrator), pointer :: prop
        double precision , allocatable:: FX(:) 

        integer :: i

        N3 = size(r0)
        Natom = N3/3
        allocate(LJA,source=LJA_)
        allocate(LJC,source=LJC_)

        allocate(Glin(N3**2), Hnew(N3, N3), GH(N3, N3), FX(N3))
        G(1:N3, 1:N3) => Glin(1:N3**2)

        kT = kT_

        allocate (lsode :: prop)

        call prop % init(N3**2, 0d0, 1d-2)
        call prop % set_dtmin(0d0)
        call prop % set_dtmax(100d0)

        prop % rtol = 1d-5
        prop % atol = 1d-4

        !G = pinv(hessian(r0))

        Glin=0d0
        Glin(1::N3+1) = 0.01d0
        call regtransrot(r0, G, 0.0d0)

        !        call frprmn(ff, fdf, r0, 1d-7,iter,F)
        do
            call prop % converge(fg, Glin, 1d-2)
            call regtransrot(r0, G, 0.00d0)
            call frprmn(ff, fdf, r0, 1d-7,iter,F)
            call regtransrot(r0, G, 0.0d0)

            print *, 'F = ', F, iter
            if (iter==1) then
                F = ff(r0)
                print *, kT, F, sqrt(sum(FX**2)/size(FX))
                exit
                if (kT>.8d0) exit
                kT = kT + 0.1
            end if
        end do

        deallocate(prop, Glin, Hnew, GH, FX)

    contains
        subroutine fg(NEQ, T, Y, YP)
            implicit none
            integer, intent(in) :: NEQ
            double precision, intent(in) :: T
            double precision, intent(in), target ::  Y(:)
            double precision, intent(out),target :: YP(:)

            double precision, pointer :: G(:,:)

            G(1:N3,1:N3) => Y(1:NEQ)

            call ForeignFreeEnergy(r0, G, F, FX, Hnew)

            call dsymm('L', 'L', N3, N3, 1d0, G, N3, Hnew, N3, 0d0, GH, N3)
            call dsymm('R', 'L', N3, N3, 1d0, G, N3, GH, N3, 0d0, yp, N3)
            yp = y - yp
            G(1:N3,1:N3) => YP(1:NEQ)
            print *, 'conv'
            call regtransrot(r0, G, 0.0d0)
        end subroutine fg
    end subroutine scp

    function ff(q)
        double precision, intent(in) :: q(:)
        double precision :: FF

        call ForeignFreeEnergy(q, G, FF)
        FF= FF/Natom
    end function ff


    subroutine fdf(q, F, FX)
        implicit none
        double precision, intent(in) :: q(:)
        double precision, intent(out) :: F
        double precision, intent(out), allocatable :: FX(:)

        if (.NOT. allocated(FX)) then
            allocate(FX(size(q)))
        end if

        call ForeignFreeEnergy(q, G, F, FX)
        F = F/Natom
        FX = FX/Natom
    end subroutine fdf


    subroutine ForeignFreeEnergy(q, G, F, UX, UXY)
        !    use omp_lib
        IMPLICIT NONE
        double precision, intent(in) :: q(:), G(:,:)
        double precision, intent(out) :: F
        double precision, intent(out), optional :: UX(:), UXY(:,:)
        INTEGER :: I1,I2,IG, J
        double precision :: AG(3,3), &
                DETA,DETAG,QZQ,U12, &
                G12(3,3),A(3,3), &
                Zq(3), Z(3,3),Q12(3)
        real*8 :: UXY0(3,3), UX0(3)

        F = 0

        if (present(UX)) UX = 0
        if (present(UXY)) UXY = 0

        do I1 = 1, N3-3, 3
            do I2 = I1+3, N3, 3
                Q12 = q(I1 : I1+2) - q(I2 : I2+2)

                G12 = 2d0 * kT * g_cross_block(G, I1, I2)

                call detminvm(G12, DETA, A)
                DETA = 1.0d0/DETA

                UX0 = 0d0; UXY0 = 0d0
                DO IG=1, size(LJC)      ! BEGIN SUMMATION OVER GAUSSIANS
                    AG = A
                    do J=1,3
                        AG(J,J) = LJA(IG) + AG(J,J)
                    end do

                    call detminvm(AG, DETAG, Z)
                    Z = - LJA(IG)**2 * Z

                    do J=1,3
                        Z(J,J) = Z(J,J) + LJA(IG)
                    end do

                    Zq = matmul(Z, Q12) ! R = -2.0*Zq
                    qZq = dot_product(Q12, Zq) 

                    if (DETA/DETAG <0) then
                        print *, 'ltzero'
                    endif
                    U12 = SQRT(DETA/DETAG) * EXP(-qZq) * LJC(IG)

                    F = F + U12

                    if (present(UX)) then
                        UX0 = UX0 - 2d0*U12*Zq
                    end if
                    if (present(UXY)) then
                        do J=1,3
                            UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                        end do
                    end if
                end do ! IG

                if (present(UX)) then
                    UX (I1 : I1+2) = UX (I1 : I1+2) + UX0
                    UX (I2 : I2+2) = UX (I2 : I2+2) - UX0
                end if

                if (present(UXY)) then
                    UXY0 = 0.5d0 * (UXY0 + transpose(UXY0))

                    UXY (I1 : I1+2, I1 : I1+2) =  UXY (I1 : I1+2, I1 : I1+2) + UXY0
                    UXY (I2 : I2+2, I2 : I2+2) =  UXY (I2 : I2+2, I2 : I2+2) + UXY0

                    UXY (I1 : I1+2, I2 : I2+2) = -UXY0
                    UXY (I2 : I2+2, I1 : I1+2) = -UXY0
                end if
            end do ! I2
        end do ! I1

        F = F - 0.5d0 * kT * logdetG(q, G) &
                - N3 * kT * log(2d0 * M_PI * exp(0.5) * kT)
    end subroutine ForeignFreeEnergy


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


    function pinv(A)
        double precision :: A(:,:)
        double precision :: pinv(size(A, 1),size(A, 1))
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


    function logdetG(q, A)
        double precision :: q(:), A(:,:)
        double precision :: logdetG
        double precision, allocatable :: AA(:,:), work(:), W(:), U(:,:)
        integer,allocatable :: iwork(:),isuppz(:)

        integer :: N, i, info, NEV

        N = size(A, 1)

        allocate(AA,source=A)
        call regtransrot(q, AA, 1d0)
        allocate(work(N * 26), W(N), iwork(N * 10), isuppz(2*N))
        call dsyevr('N', 'A', 'L', N, AA, N, 0d0, 0d0, 1, 6, &
                DLAMCH('Safe minimum'), NEV, W, 0, N, isuppz, &
                work, size(work), iwork, size(iwork), info)

        !print '("EV:", 12F12.6)',W(1:6),W(N-5:N)
        logdetG = sum(log(abs(W)))

        deallocate(AA,W,work,iwork,isuppz)
    end function logdetG


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

    function hessian(q) result(UXY)
        !    use omp_lib
        IMPLICIT NONE
        double precision, intent(in) :: q(:)
        double precision :: UXY(size(q), size(q))

        double precision  :: U12, Q12(3), UXY0(3,3)
        integer :: J,I1,I2,IG

        UXY = 0d0
        do I1=1,size(q),3
            DO I2=I1+3,size(q),3
                Q12 = q(I1 : I1+2) - q(I2 : I2+2)

                UXY0 = 0d0
                DO IG=1,size(LJA)
                    U12 = EXP(-LJA(IG) * sum(Q12**2) ) * LJC(IG)

                    do J=1,3
                        UXY0(:,J) = UXY0(:,J) + 4d0*U12*LJA(IG)**2 * Q12(j) * Q12
                    end do

                    do J=1,3
                        UXY0(J,J) = UXY0(J,J) - 2d0*U12 * LJA(IG)
                    end do
                end do
                ! guarantee symmetry
                UXY0 = 0.5d0 * (UXY0 + transpose(UXY0))

                UXY(I1 : I1+2, I1 : I1+2) = UXY(I1 : I1+2, I1 : I1+2) + UXY0
                UXY(I2 : I2+2, I2 : I2+2) = UXY(I2 : I2+2, I2 : I2+2) + UXY0

                UXY(I1 : I1+2, I2 : I2+2) = -UXY0
                UXY(I2 : I2+2, I1 : I1+2) = -UXY0
            end do ! I2
        end do ! I1
    END function hessian

    function transrot_subspace(q) result(U)
        implicit none
        double precision, intent(in) :: q(:)
        double precision :: U(size(q),6)

        integer :: i, j, N

        N = size(q)

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
    end function transrot_subspace

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

    subroutine slash_g(q, G)
        double precision :: q(:), G(:,:)

        double precision, allocatable :: U(:,:), GU(:,:)
        integer :: i, j, N

        N = size(G, 1)

        print *, 'slash!!'

        allocate(U(N,6), GU(N,6))

        U = transrot_subspace(q)

        call dsymm('L', 'L', N, 3, 1d0, G, N, U, N, 0d0, GU, N)
        call dgemm('N', 'T', N, N, 3, -1d0, U, N, GU, N, 1d0, G, N)
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
end module scp_mod
