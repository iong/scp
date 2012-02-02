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
        integer :: N3, info, i

        double precision, pointer :: G(:,:)

        self % NG =  9 * self%Natom**2
        self % NEQ = 3 * self % Natom + self % NG

        N3 = 3 * self % Natom
        if (.NOT. allocated(self % y)) then
            allocate(self % y ( self % NEQ ), self%UXY(N3, N3))
        end if

        allocate (lsode :: self%prop)

        call self % prop % init(self%NEQ, 0d0, self%dt0)
        call self % prop % set_dtmin(0d0*self%dt_min)
        !call self % prop % set_dtmax(1d-2)

        self % prop % rtol = 1d-5
        self % prop % atol (1 : 3 * self % Natom) = 1d-3!self % q_atol
        self % prop % atol (3 * self % Natom + 1 : ) = 1d-5!self % g_atol

        self % y(1 : N3) = reshape(Q0, (/ N3 /) )

        G(1:N3,1:N3) => self%y(N3+1:)
        call self % hessian( self % y(1 :N3), G)
        call dpotrf('L', N3, G, N3, info)
        call dpotri('L', N3, G, N3, info)

        !G = 0d0
        do i=1,N3
            G(i,i:) = G(i:,i)
            !G(i,i) = 1d0
        end do
    end subroutine

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
            INTEGER :: J,I1,I2,IG, I1_3, I2_3, skip, N3, info
            REAL*8 AG(3,3), &
                    DETA,DETAG,QZQ,U12, &
                    G12(3,3),A(3,3), &
                    Zq(3), Z(3,3),Q12(3)
            real*8 :: UXY0(3,3), UX0(3), U, TRUXXG, k1

            double precision, pointer :: UX(:), GPP(:), G(:,:), GP(:,:)

           
            N3 = 3*self%Natom
            UX => yp(1:N3)
            G(1:N3,1:N3) => y(N3+1:)
            GP(1:N3,1:N3) => yp(N3+1:)

            self % U = 0; UX = 0; self%UXY = 0;
            do I1 = 1, self%Natom - 1
                I1_3 = 3*(I1-1) + 1
                do I2 = I1+1, self%Natom
                    I2_3 = 3*(I2-1) + 1
                    Q12 = y(I1_3:I1_3+2) - y(I2_3:I2_3+2)
                    Q12 = min_image(Q12, self%bl)
                    G12=  G (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) &
                            + G (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) &
                            - G (I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) &
                            - G (I1_3 : I1_3 + 2, I2_3 : I2_3 + 2)
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
                    UX (I1_3 : I1_3 + 2) = UX (I1_3 : I1_3 + 2) + UX0
                    UX (I2_3 : I2_3 + 2) = UX (I2_3 : I2_3 + 2) - UX0

                    self%UXY (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) = &
                            self%UXY (I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) + UXY0

                    self%UXY (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) = &
                            self%UXY (I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) + UXY0

                    !if (sum(q12**2) <= rfullmatsq) then
                    self%UXY (I1_3 : I1_3 + 2, I2_3 : I2_3 + 2) = -UXY0
                    self%UXY (I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) = -UXY0
                    !end if
                end do ! I2
            end do ! I1

            call dsymm('L', 'L', N3, N3, 1d0, G, N3, self%UXY, N3, 0d0, GP, N3)
            do i1=1,N3
                GP(i1,i1) = GP(i1,i1) - 1d0
            end do
            GP = 0.5d0  *  ( GP + transpose(GP) )

            self % qconv = sqrt(sum(UX**2)/size(UX))
            self % gconv = sqrt(sum(GP**2)/size(GP))

            UX = -self%cq * UX
            GP = -0.5d0* self%cg * self%kT * GP
        END SUBROUTINE RHS
    end subroutine PROPAGATE

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

        double precision, pointer :: G(:,:)

        integer :: info, j

        G(1:3 * self%Natom, 1:3 * self%Natom) => self%y(3 * self%Natom + 1:)
        call dpotrf('U', 3 * self%Natom, G, 3 * self%Natom, info)

        logdet=0.0
        DO j=1, 3 * self%Natom
            logdet = logdet + LOG(ABS( G(j,j) ))
        ENDDO
        logdet = 2d0* logdet
    end function logdet

end module vgwfm_mod
