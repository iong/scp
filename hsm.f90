program scp
    use utils, only: M_PI
    use xyz, only: load_xyz
    implicit none


    character(256) :: fname

    integer :: i, Natom

    double precision, allocatable :: r0(:,:), H(:,:)


    call get_command_argument(1, fname)
    call load_xyz(fname, r0)
       
    Natom = size(r0, 2)

    allocate(H(3*Natom, 3*Natom))

    H = hessian(r0)

    fname = fname(1:index(fname,'.xyz')-1) // '_spectr.dat'
    call printev(H, fname)
 
    deallocate(r0, H)
contains
    subroutine printev(A, fname)
        double precision :: A(:,:)
        character(*) :: fname
        double precision, allocatable :: work(:), W(:), U(:,:)
        integer,allocatable :: iwork(:),isuppz(:)

        integer :: N, info, NEV
        interface
            double precision function dlamch(s)
                character(1) :: s
            end function dlamch
        end interface

        N = size(A, 1)

        allocate(U(N,N))
        allocate(work(N * 26), W(N), iwork(N * 10), isuppz(2*N))
        call dsyevr('V', 'A', 'L', N, A, N, 0d0, 0d0, 1, 6, &
                DLAMCH('Safe minimum'), NEV, W, U, N, isuppz, &
                work, size(work), iwork, size(iwork), info)
    
        open(35, file=trim(fname))
        write(35,'(F18.10)') W
        close(35)

        deallocate(W,work,iwork,isuppz, U)
    end subroutine printev

    function hessian(q) result(UXY)
        !    use omp_lib
        IMPLICIT NONE
        double precision, intent(in) :: q(:,:)
        double precision :: UXY(size(q), size(q))

        double precision  :: U12, Q12(3), UXY0(3,3)
        double precision :: LJA(3) = (/ 6.65, 0.79, 2.6 /), &
                LJC(3) = (/ 1840d0, -1.48d0, -23.2d0 /)
        
        integer :: J,I1,I2,IG

        print *, size(q)
        UXY = 0d0
        do I1=1,3*size(q,2),3
            DO I2=I1+3,3*size(q,2),3
                Q12 = q(:,(I1+2)/3) - q(:,(I2+2)/3)

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
end program scp
