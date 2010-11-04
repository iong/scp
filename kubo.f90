subroutine kubo(q0, v0, beta, nsteps, xk, vk)
    use spine
    use utils
    use vgw
    implicit none
    real*8, intent(in) :: q0(:,:), v0(:,:), beta
    integer, intent(in) :: nsteps
    real*8, intent(out) :: xk(:,:), vk(:,:)
    real*8, dimension(3, Natom) :: ql, qr, dxk, dvk
    real*8, dimension(3, 3, Natom) :: gl, gr, Qnkl, Qnkr
    real*8, dimension(1+21*Natom, 0:nsteps) :: yl
    real*8 :: tau, dtau, Ueff, dqlr(3), invglr(3,3), logdet, qGq, &
            vnkl(3), vnkr(3), detj, Ueffbeta, w, qlqr, qlqrint
    integer :: i, j
 
    if (mod(nsteps,2) /= 0) then
        write (*,*) 'nsteps must be even!!'
        stop
    end if

    tau = 0.0d0
    dtau = beta / nsteps
    call init_gaussians(q0, 0.0d0, yl(:,0))
    do i=1,nsteps
        yl(:,i) = yl(:,i-1)
        call vgw1(q0, Ueff, tau+2.0*dtau, tau, yl(:,i))
        if (i == nsteps/2) then
            Ueffbeta  = Ueff
        end if
        tau = tau+2.0*dtau
    end do

    xk = 0.0d0
    vk = 0.0d0
    qlqrint = 0.0d0
    do i=0,nsteps/2-1
        call unpack_q(yl(:,i), ql)
        call unpack_g(yl(:,i), gl)
        call unpack_qnk(yl(:,i), Qnkl)
        call unpack_q(yl(:,nsteps-i), qr)
        call unpack_g(yl(:,nsteps-i), gr)
        call unpack_qnk(yl(:,nsteps-i), Qnkr)
        if (i==0) then
            w = 1.0d0
        else
            w = 2.0d0
        end if

        qGq = 0.0d0
        logdet = 0.0d0
        do j=1,Natom
            dqlr = ql(:,j) - qr(:,j)
            call detminvm(gl(:,:,j) + gr(:,:,j), detj, invglr)
            qGq = qGq + dot_product(dqlr, matmul(invglr, dqlr))
            logdet = logdet + log(detj)
            dxk(:,j) = matmul(invglr, matmul(gr(:,:,j), ql(:,j)) + matmul(gl(:,:,j), qr(:,j)))
            vnkl = matmul(Qnkl(:,:,j), v0(:,j))
            vnkr = matmul(Qnkr(:,:,j), v0(:,j))
            dvk(:,j) = matmul(invglr, matmul(gr(:,:,j), vnkl) + matmul(gl(:,:,j), vnkr))
        end do
        qlqr = exp(yl(1, i) + yl(1, nsteps-i) - 0.5*(qGq + logdet))
        qlqrint = qlqrint + w*qlqr
        xk = xk + dxk*w
        vk = vk + dvk*w
    end do
    write (*,*) beta, (qlqrint * exp(beta*Ueffbeta + 1.5d0*Natom*log(2.0d0)) + 1.0d0)/nsteps
    !write (*,*) 'qlqrint =',qlqrint/nsteps * exp(beta*Ueffbeta + 1.5d0*Natom*log(2.0d0))

    call unpack_q(yl(:,nsteps/2), ql)
    call unpack_qnk(yl(:,nsteps/2), Qnkl)
    do j=1,Natom
        dvk(:,j) = matmul(Qnkl(:,:,j), v0(:,j))
    end do
    
    xk = (xk + ql) / nsteps
    vk = (vk + dvk) / nsteps
end subroutine
