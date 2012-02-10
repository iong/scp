program gaussmin
    use nr, only: frprmn
    use utils, only : M_PI, detminvm
    use xyz
    double precision, allocatable :: r0(:), H(:,:), Hnew(:,:)
    integer :: Natom, i, iter=1000
    character(256) :: coords, waste
    double precision :: U0, kT
    double precision, parameter :: LJA(3) = (/ 6.65, 0.79, 2.6 /), &
          LJC(3) = (/ 1840d0, -1.48d0, -23.2d0 /)

    call get_command_argument(1, coords)

    open(33, file=trim(coords))
    read (33, *) Natom

    allocate(r0(3*Natom), H(3*Natom, 3*Natom), Hnew(3*Natom, 3*Natom))

    read(33,*)
    read(33, *) (waste, r0(3*i-2:3*i), i=1,Natom)
    close(33)

    if (command_argument_count() > 1) then
        call get_command_argument(2, waste)
        read(waste, *) kT
    end if

    call frprmn(Utot, UtotX, r0, 1d-10,iter,U0)
    print *, Natom, U0, iter

    H = hessian(r0)
    call frprmn(ff, fdf, r0, 1d-10,iter,U0)
    print *, Natom, U0, iter

    coords(len_trim(coords) - 3:) = '_gauss.xyz'
    write(waste, *) U0
    call dump_xyz(reshape(r0, (/3, Natom/)), coords, waste)
    
contains
    function Utot(q)
        implicit none
        double precision, intent(in) :: q(:)
        double precision :: Utot

        integer :: i1, i2
        double precision :: q12(3)

        Utot = 0d0
        do I1=1,size(q),3
            DO I2=I1+3,size(q),3
                Q12 =q(I1 : I1+2) - q(I2 : I2+2)

                Utot = Utot  + sum(EXP(-LJA * sum(Q12**2) ) * LJC)
            end do ! I2
        end do ! I1
    end function Utot

    subroutine UtotX(q, U, UX)
        implicit none
        double precision, intent(in) :: q(:)
        double precision, intent(out) :: U
        double precision, allocatable, intent(out) :: UX(:)

        integer :: i1, i2
        double precision :: q12(3), UX0(3)

        if (.NOT. allocated(UX)) allocate(UX(size(q)))

        U = 0d0
        UX = 0d0
        do I1=1,size(q),3
            DO I2=I1+3,size(q),3
                Q12 =q(I1 : I1+2) - q(I2 : I2+2)
                U = U  + sum(EXP(-LJA * sum(Q12**2) ) * LJC)
                ! there are to -1 here which cancel each other
                UX0 = -q12*(2.0*sum( LJA * LJC * EXP(-sum(q12**2) * LJA) ) )

                UX(I1:I1+2) = UX(I1:I1+2) + UX0
                UX(I2:I2+2) = UX(I2:I2+2) - UX0
            end do ! I2
        end do ! I1
    end subroutine UtotX

    function ff(q)
        double precision, intent(in) :: q(:)
        double precision :: FF

        double precision, allocatable :: FX(:)

        call fdf(q, ff, FX)
        deallocate(FX)
    end function ff

    subroutine fdf(q, F, FX)
        implicit none
        double precision, intent(in) :: q(:)
        double precision, intent(out) :: F
        double precision, intent(out), allocatable :: FX(:)

        double precision :: dHnorm

        if (.NOT. allocated(FX)) then
            allocate(FX(size(q)))
        end if

        do
            call ForeignFreeEnergy(q, H, F, FX, Hnew)
             dHnorm = sqrt(sum((Hnew - H)**2)/size(H))
             print *, 'dHnorm = ', dHnorm
             if (dHnorm < 1d-2) then
                 exit
             end if
             H = Hnew
         end do
     end subroutine fdf

    subroutine ForeignFreeEnergy(q, UXYold, F, UX, UXY)
        !    use omp_lib
        IMPLICIT NONE
        double precision, intent(in) :: q(:), UXYold(:,:)
        double precision, intent(out) :: F, UX(:), UXY(:,:)
        INTEGER :: I1,I2,IG, J
        double precision :: AG(3,3), &
              DETA,DETAG,QZQ,U12, &
              G12(3,3),A(3,3), &
              Zq(3), Z(3,3),Q12(3)
        real*8 :: UXY0(3,3), UX0(3), logdetG

        double precision, allocatable :: G(:,:)

        allocate(G(3*Natom,3*Natom))
        call pinv(UXYold, G, logdetG)

        UXY = 0
        UX = 0
        F = 0
        do I1 = 1, 3*(Natom-1),3
            do I2 = I1+3, 3*Natom, 3
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
  
                    UX0 = UX0 - 2d0*U12*Zq
                    do J=1,3
                        UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                    end do
                end do ! IG

                UX (I1 : I1+2) = UX (I1 : I1+2) + UX0
                UX (I2 : I2+2) = UX (I2 : I2+2) - UX0

                UXY0 = 0.5d0 * (UXY0 + transpose(UXY0))

                UXY (I1 : I1+2, I1 : I1+2) =  UXY (I1 : I1+2, I1 : I1+2) + UXY0
                UXY (I2 : I2+2, I2 : I2+2) =  UXY (I2 : I2+2, I2 : I2+2) + UXY0

                UXY (I1 : I1+2, I2 : I2+2) = -UXY0
                UXY (I2 : I2+2, I1 : I1+2) = -UXY0
            end do ! I2
        end do ! I1

        F = F - 0.5d0 * kT * logdetG &
                - 3 * Natom * kT * log(2d0 * M_PI * exp(0.5) * kT)
        deallocate(G)
    END subroutine ForeignFreeEnergy

    subroutine fill_lower(y, A)
        implicit none
        double precision, intent(in) :: y(:)
        double precision, intent(out) :: A(:,:)

        integer :: j, i, N

        N = size(A, 1)
        j = 1
        A=0d0
        do i = 1, N
            A(i:,i) = y(j : j + N - i)
            j = j + N - i + 1
        end do
    end subroutine fill_lower

    subroutine pack_lower(A, y)
        implicit none
        double precision, intent(in) :: A(:,:)
        double precision, intent(out) :: y(:)

        integer :: j, i, N

        N = size(A, 1)
        j = 1
        do i = 1, N
            y(j : j + N - i) = A(i:,i)
            j = j + N - i + 1
        end do
    end subroutine pack_lower

    pure function g_cross_block(G, I, J) result(x)
        double precision, intent(in) :: G(:,:)
        integer, intent(in) :: I, J
        double precision ::x(3,3)
        x = G(I : I+2, I : I+2) + G(J : J+2, J : J+2)
        x(1,2) = x(2,1)
        x(1,3) = x(3,1)
        x(2,3) = x(3,2)
        if (I>J) then
            x = x - G(I : I+2, J : J+2) - transpose( G(I : I+2, J : J+2) )
        else
            x = x - G(J : J+2, I : I+2) - transpose( G(J : J+2, I : I+2) )
        end if
    end function g_cross_block

    function logdet(A)
        double precision, intent(in) :: A(:,:)

        integer :: info, j, N

        N = size(A,1)
        call dpotrf('L', N, A, N, info)

        logdet=0.0
        DO j=1, N
            logdet = logdet + LOG(ABS( A(j,j) ))
        ENDDO
        logdet = 2d0* logdet
    end function logdet

    subroutine pinv(A, Ainv, logdetAinv)
        double precision :: A(:,:)
        double precision, intent(out) :: Ainv(:,:), logdetAinv
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
        call dgemm('N', 'T', N, N, N, 1d0, U, N, AA, N, 0d0, Ainv, N)

        logdetAinv = -sum(log(abs(W(7:))))

        deallocate(AA,W,U,work,iwork,isuppz)
    end subroutine pinv

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
end program gaussmin
