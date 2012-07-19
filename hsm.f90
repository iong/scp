program scp
    use utils, only: M_PI
    use xyz, only: load_xyz
    implicit none


    character(256) :: fname

    integer :: i, Natom

    double precision :: U0
    double precision, allocatable :: r0(:,:), H(:,:), W(:)


    call get_command_argument(1, fname)
    call load_xyz(fname, r0)
       
    Natom = size(r0, 2)

    allocate(H(3*Natom, 3*Natom), W(3*Natom))

    call ghessian(r0, U0, H)
    W = getev(H)
    deallocate(r0, H)

    fname = fname(1:index(fname,'.xyz')-1) // '_spectr.dat'
    open(35, file=trim(fname))
    write(35,'(F18.10)') W
    close(35)
    print *, Natom, U0/Natom, sum(sqrt(W(7:)))/(3*Natom - 6), &
        exp(0.5d0*sum(log(W(7:)))/(3*Natom - 6))
contains
    function getev(A) result(W)
        double precision, intent(in) :: A(:,:)
        double precision :: W(size(A, 1))
        double precision, allocatable :: work(:), U(:,:)
        integer,allocatable :: iwork(:),isuppz(:)

        integer :: N, info, NEV
        interface
            double precision function dlamch(s)
                character(1) :: s
            end function dlamch
        end interface

        N = size(A, 1)

        allocate(U(N,N))
        allocate(work(N * 26), iwork(N * 10), isuppz(2*N))
        call dsyevr('V', 'A', 'L', N, A, N, 0d0, 0d0, 1, 6, &
                DLAMCH('Safe minimum'), NEV, W, U, N, isuppz, &
                work, size(work), iwork, size(iwork), info)
    


        deallocate(work,iwork,isuppz, U)
    end function getev

    subroutine ghessian(q, Utot, UXY)
        IMPLICIT NONE
        double precision, intent(in) :: q(:,:)
        double precision, intent(out) :: Utot, UXY(size(q), size(q))

        double precision  :: U12, Q12(3), UXY0(3,3)
!        double precision :: LJA(3) = (/ 6.65, 0.79, 2.6 /), &
!                LJC(3) = (/ 1840d0, -1.48d0, -23.2d0 /)
        double precision :: LJA(4) = (/13.75241581543,  7.64136762716565, &
                1.66106514430194,  0.46671857809601/)
        double precision :: LJC(4) = (/128891.569912461, 3513.506510626, &
                -8.47940219671101, -0.343373800123181/)
        
        integer :: J,I1,I2,IG

        print *, size(q)
        Utot = 0d0
        UXY = 0d0
        do I1=1,3*size(q,2),3
            DO I2=I1+3,3*size(q,2),3
                Q12 = q(:,(I1+2)/3) - q(:,(I2+2)/3)

                UXY0 = 0d0
                DO IG=1,size(LJA)
                    U12 = EXP(-LJA(IG) * sum(Q12**2) ) * LJC(IG)

                    Utot = Utot + U12

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
    END subroutine ghessian

    subroutine hessian(q, Utot, UXY)
        !    use omp_lib
        IMPLICIT NONE
        double precision, intent(in) :: q(:,:)
        double precision,intent(out) :: Utot, UXY(size(q), size(q))

        double precision  :: U12, Q12(3), UXY0(3,3), qsq, y8, y6, A, B        
        integer :: J,I1,I2

        Utot = 0d0
        UXY = 0d0
        do I1=1,3*size(q,2),3
            DO I2=I1+3,3*size(q,2),3
                Q12 = q(:,(I1+2)/3) - q(:,(I2+2)/3)

                qsq = sum(q12**2)
                y8 = 1d0/qsq**4

                y6 = y8 * qsq

                Utot = Utot + 4d0*y6*(y6 - 1d0)

                A = (672d0*y8 - 192d0/qsq)*y8;
                do J=1,3
                    UXY0(:,J) = (A*Q12(J))*Q12
                end do

                B = 24d0*y8 - 48d0*y8*y8*qsq
                do J=1,3
                    UXY0(J,J) = UXY0(J,J) + B
                end do
                
                UXY(I1 : I1+2, I1 : I1+2) = UXY(I1 : I1+2, I1 : I1+2) + UXY0
                UXY(I2 : I2+2, I2 : I2+2) = UXY(I2 : I2+2, I2 : I2+2) + UXY0

                UXY(I1 : I1+2, I2 : I2+2) = -UXY0
                UXY(I2 : I2+2, I1 : I1+2) = -UXY0
            end do ! I2
        end do ! I1
    END subroutine hessian
end program scp
