SUBROUTINE RHSSspFM(NEQ, T, Y, YP)!, RPAR, IPAR)
    use sparse
    IMPLICIT NONE
    integer, intent(in) :: NEQ!, IPAR(:)
    double precision, intent(in) :: T!, RPAR(:)
    double precision, intent(in), target :: Y(:)
    double precision, intent(out), target :: YP(:)

    INTEGER :: J,I1,I2,IG, p
    double precision :: AG(3,3), DETA,DETAG,QZQ,U12, Gb12(3,3),A(3,3), &
            Zq(3), Z(3,3),Q12(3), UXY0(3,3), UX0(3)

    type(csr) :: GUG
    
    integer :: Gbandpos, UXYbandpos

    if (y(3*Natom+1)==0d0) then
        call rhss_zero_time(y, yp)
        return
    end if

    G%x => y(3*Natom+1 : 3*Natom + G%nnz)

    call G%get_3x3_diag(Gbdiag)

    U = 0; UX = 0;  UXY%x=0; UXYbdiag=0;
!$omp parallel default(private) shared(G, UXY, y, Gbdiag, Natom, NGAUSS, LJA, LJC, LAMBDA) reduction(+:U,UX,UXYbdiag)
    
    allocate( Gband(3,G%nnz_row_max), Gband_ja(G%nnz_row_max),&
          UXYband(3, UXY%nnz_row_max))
 
!$omp do schedule(dynamic)
    do I1=1,3*(Natom-1), 3
        Gbandpos = 1
        Gband_ja = 0
        if (G%iia(I1) + 3 < G%ia(I1+1)) then
            call G%get_3rows_ur(I1, Gband)
            Gband_ja(:G%ia(I1+1) - G%iia(I1) - 3) = G%ja(G%iia(I1) + 3 : G%ia(I1 + 1) - 1)
        end if

        UXYband = 0d0

        do p=UXY%iia(I1) + 3, UXY%ia(I1+1)-1, 3 
            I2 = UXY%ja(p)

            Q12 = y(I1 : I1 + 2) - y(I2 : I2 + 2)
            Q12 = min_image(Q12, bl)

            Gb12=Gbdiag(:,I1 : I1+2) + Gbdiag(:,I2 : I2+2)
            if ( Gbandpos < G%nnz_row_max .and. Gband_ja(Gbandpos) == I2) then

                Gb12 = Gb12 - Gband(:,Gbandpos : Gbandpos + 2) &
                        - transpose(Gband(:,Gbandpos : Gbandpos + 2))

                Gbandpos = Gbandpos + 3
            end if

            Gb12 = Gb12 * LAMBDA

            call detminvm(Gb12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXY0 = 0d0
            DO IG=1,NGAUSS ! BEGIN SUMMATION OVER GAUSSIANS
                AG = A
                do J=1,3
                    AG(J,J)=LJA(IG)+AG(J,J)
                end do

                call detminvm(AG, DETAG, Z)
                Z = - LJA(IG)**2 * Z

                do J=1,3
                    Z(J,J) = Z(J,J) + LJA(IG)
                end do

                Zq = matmul(Z, Q12) ! R = -2.0*Zq
                qZq = dot_product(Q12, Zq) 

                if (DETA*DETAG <= 0.0d0 ) then
                    write (69, '(F16.8)'), Y
                    write (*,*) DETA, DETAG
                    stop
                end if
                U12 = SQRT(DETA/DETAG)*EXP(-qZq)*LJC(IG)
                U = U + U12

                UX0 = UX0 - 2d0*U12*Zq
                do J=1,3
                    UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                end do
            end do ! IG

            UX(I1 : I1 + 2) = UX(I1 : I1 + 2) + UX0
            UX(I2 : I2 + 2) = UX(I2 : I2 + 2) - UX0

            UXYbdiag(:,I1 : I1 + 2) = UXYbdiag(:,I1 : I1 + 2) + UXY0
            UXYbdiag(:,I2 : I2 + 2) = UXYbdiag(:,I2 : I2 + 2) + UXY0

            UXYbandpos = P - UXY%iia(I1) + 1
            UXYband(:,UXYbandpos : UXYbandpos + 2) = -UXY0

        end do ! I2
        call UXY%write_3rows_urnd(I1, UXYband(:,4:))
      
    end do ! I1
!$omp end do
    deallocate (Gband_ja, Gband, UXYband)
!$omp end parallel
    call UXY%set_3x3_diag(UXYbdiag)
    call UXY%mirror_uplo()

    call G % gemv(UX, yp)
    yp(1:3*Natom) = -yp(1:3*Natom)
 
    call G % multiply_restricted(UXY, GU)

    GUG = G
    GUG%x => yp(3*Natom+1 : 3*Natom + G%nnz)
    call GU % multiply_restricted(G, GUG)
    GUG%x = -GUG%x
    GUG%x (GUG%iia) = GUG%x (GUG%iia) + invmass
    
    call GUG % force_symmetry()

    yp(NEQ) = -(0.25d0 * LAMBDA * GU%trace() + U)/real(Natom)

! full matrix debug
!!$    call dsymm('L', 'L', 3*Natom, 1, -1d0, G%x, 3*Natom, UX, 3*Natom, &
!!$        0d0, GP, 3*Natom)
!!$
!!$    !write (*,'(F16.8$)') maxval(abs(GP(1:3*Natom) - yp(1:3*Natom)))
!!$
!!$    call dsymm('L', 'L', 3*Natom, 3*Natom, 1d0, G%x, 3*Natom, UXY%x, 3*Natom, &
!!$        0d0, GU%x, 3*Natom)
!!$    call dsymm('R', 'L', 3*Natom, 3*Natom, -1d0, G%x, 3*Natom, GU%x, 3*Natom, &
!!$        0d0, GP, 3*Natom)
!!$
!!$    GP(::3*Natom+1) = GP(::3*Natom+1) + invmass
!!$
!!$    !print *, maxval(abs(GP - GUG%x))
!!$
!!$    GUG%x = GP
!!$    yp(NEQ) = -(0.25d0 * GU%trace() + U)/real(Natom)

    
    nullify(G%x)
end subroutine RHSSspFM

subroutine rhss_zero_time(y, yp)
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), rsq
    integer :: i, j, Jp

    yp = 0
    yp(3*Natom + G%iia) = invmass

    U=0d0

    DO I=1,UXY%nrows-3,3
        DO Jp=UXY%iia(I) + 3, UXY%ia(I+1) - 1, 3
            J = UXY%ja(Jp)
            qij = y(I : I+2) - y(J : J+2)
            rsq = sum(min_image(qij, BL)**2)
            U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO

    yp(3*Natom + G%nnz + 1) = -U/real(Natom)
end subroutine rhss_zero_time
