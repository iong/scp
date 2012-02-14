program scpmin
    use nr, only: frprmn
    use scp_mod
    use xyz
    implicit none
    double precision, allocatable :: r0(:)
    integer :: Natom, N3, i, iter=1000
    character(256) :: coords, waste
    double precision :: U0, kT
    double precision, parameter :: LJA(3) = (/ 6.65, 0.79, 2.6 /), &
          LJC(3) = (/ 1840d0, -1.48d0, -23.2d0 /)

    call get_command_argument(1, coords)

    open(33, file=trim(coords))
    read (33, *) Natom

    N3 = 3*Natom
    allocate(r0(N3))

    read(33,*)
    read(33, *) (waste, r0(3*i-2:3*i), i=1,Natom)
    close(33)

    if (command_argument_count() > 1) then
        call get_command_argument(2, waste)
        read(waste, *) kT
    end if

    call frprmn(Utot, UtotX, r0, 1d-10,iter,U0)
    print *, Natom, U0, iter

    call scp(r0, LJC, LJA, kT, U0, iter)
   
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
end program scpmin
