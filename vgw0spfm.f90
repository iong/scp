SUBROUTINE vgw0spfm_init_prop(Q0, BL_)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), BL_

    double precision, allocatable :: q(:,:)
    integer :: nnz(2)

    type(csr_list) :: matrix_list(2)

    Natom = size(Q0, 2)
    BL = BL_ / sigma0
    allocate(q, source=q0/sigma0)

    nnz = count_pairs(q, (/ Vcutoff, Gcutoff /) )
    nnz = 9*(2*nnz + Natom)
    call UXY%init(3*Natom, 3*Natom, nnz(1))
    call G%init(3*Natom, 3*Natom, nnz(2))

    matrix_list(1)%p => UXY
    matrix_list(2)%p => G
    call init_sparse_patterns(Q, (/ Vcutoff, Gcutoff /), matrix_list )

    GU = G%multiply(UXY)
    call GU % sort()
    call GU%find_diagonal()

    if (debug >= info) then
        print '("G  sparsity = ",F7.3,"%")', 100d0 * G%nnz / ( G%nrows * G%ncols)
        print '("U  sparsity = ",F7.3,"%")', 100d0 * UXY%nnz / ( UXY%nrows * UXY%ncols)
        print '("GU sparsity = ",F7.3,"%")', 100d0 * GU%nnz / ( GU%nrows * GU%ncols)
    end if

    deallocate(G%x)

    NEQ = 3*Natom + G%nnz + 1

    allocate(Gbdiag(3,3*Natom), UXYbdiag(3,3*Natom), UX (3*Natom), Y(NEQ))
    
    allocate(lsode :: prop)
    call prop % init(NEQ, 0d0, dt0)
    call prop % set_dtmin(dtmin)
    call prop % set_dtmax(dtmax)

    prop % RTOL=1d-5
    prop % ATOL(1:3*Natom) = vgw_atol(1)
    prop % ATOL(3*Natom+1:3*Natom+G%nnz)=vgw_atol(2)
    prop % ATOL(NEQ) = vgw_atol(3)

    y=0d0
    y(1:3*Natom) = reshape(Q, (/ 3*Natom /) )
	
	deallocate(q)
end subroutine

SUBROUTINE vgw0spfm(Q0, BL_, beta, Ueff, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    double precision, intent(out), optional :: rt

    real*8 :: t0, t1, t2, m0, taustop

    !m0 = minval(mass)
    m0 = mass
    LAMBDA = 1d0/(sigma0 * sqrt(m0 * epsilon0))
    invmass = m0/mass
    
    call vgw0spfm_init_prop(q0, bl_)
 
    call cpu_time(t0)

    taustop = 0.5d0 * beta * sigma0 * LAMBDA
    call prop % advance(RHSSspFM, y, taustop)

    call cpu_time(t1)

    if (debug >= info) then
        write (*,*) prop % nsteps, 'steps,', prop % ncalls, ' RHSS calls'
    end if

    Ueff = -logrho()/beta
    nullify(G%x)

    call cpu_time(t2)

    if (present(rt)) then
        rt = (t2 - t0) / real(prop % ncalls)
     end if
END SUBROUTINE

double precision function logrho()
    double precision :: gama
    gama = Y(NEQ) * real(Natom) * epsilon0 / LAMBDA
    G%x => y(3*Natom+1 : 3*Natom + G%nnz)
    logrho = 2.0*gama - 0.5*G%logdet() - 1.5d0 * Natom * log(LAMBDA*sigma0**2) &
            - 1.5*Natom*log(4.0*M_PI)
end function logrho


SUBROUTINE vgw0spfmgs(Q0, BL_, beta_LAMBDA, Havg, UXabs, betafinish)
    use iso_fortran_env
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta_LAMBDA, BL_
    double precision, intent(out) :: Havg(2), UXabs(2), betafinish(2)
    real*8 :: lrho(4), dT, t0, t1, T, UXabs0, m0, taustop

    !m0 = minval(mass)
    m0 = mass
    LAMBDA = 1d0/(sigma0 * sqrt(m0 * epsilon0))
    invmass = m0/mass
    
    call vgw0spfm_init_prop(q0/sigma0, bl_/sigma0)

    call cpu_time(t0)

    T = beta_LAMBDA * epsilon0 !* LAMBDA
    dT = 0.1d0
    call prop % advance(RHSSspFM, y, 0.5*T)
    write(ERROR_UNIT,*) 'T =',T,' <UX> =', sqrt(sum(UX**2)/(3*Natom))

    lrho(1) = logrho()

    T = T + dT
    call prop % advance(RHSSspFM, y, 0.5*T)
    lrho(2) = logrho()

    betafinish(1) = T
            
    T = T + dT*10
    call prop % advance(RHSSspFM, y, 0.5*T)
    write(ERROR_UNIT,*) 'T =',T,' <UX> =', sqrt(sum(UX**2)/(3*Natom))

    lrho(3) = logrho()

    T = T + dT
    call prop % advance(RHSSspFM, y, 0.5*T)
    lrho(4) = logrho()

    betafinish(2) = T

    if (debug >= info) then
        write (*,*) prop % nsteps, 'steps,', prop % ncalls, ' RHSS calls'
     end if

    Havg = -(lrho(2::2) - lrho(1::2)) / (dT/(epsilon0 * LAMBDA))

    call cpu_time(t1)

!!$    if (present(rt)) then
!!$        rt = (stop_time - start_time) / real(prop % ncalls)
!!$     end if
END SUBROUTINE

function count_pairs(Q, srad) result(npairs)
  implicit none
  real*8, intent(in) :: Q(:,:), srad(:)
  integer ::npairs(size(srad))
  integer :: N,I,J, NN
  real*8 :: rsq,sradsq(size(srad)),qij(3)

  N = size(Q, 2)
  sradsq=srad**2

  npairs = 0
  do I=1,N-1
      do J=I+1,N
          qij=Q(:,I)-Q(:,J)
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
subroutine init_sparse_patterns(Q0, srad, ML)
    double precision, intent(in) :: Q0(:,:), srad(:)
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
    do ib=1,Natom
        i = 3*ib - 2

        p0 = p

        do jb=1,Natom
            j = 3*jb - 2

            qij=Q0(:,ib)-Q0(:,jb)
            rsq = sum(min_image(qij, BL)**2)
            
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
        ML(k)%p % ia(3*Natom + 1) = p(k)

        if (ML(k)%p % nnz /= p(k) - 1) then
            write (*,*) 'nnz error!', k, ML(k)%p%nnz, p(k)-1
            stop
        end if
    end do
end subroutine

subroutine JAC()
end subroutine

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
end  subroutine
