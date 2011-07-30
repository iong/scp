subroutine interaction_lists(Q,N,rc,bl, QRC)
      implicit none
      integer N,I,J,K
      real*8 Q(3,N),rsq,rc,rc2,bl,bl2,qij
      logical QRC(N*N)

      rc2=rc**2
      bl2=bl/2

      do I=1,N-1
        do J=I+1,N
          rsq=0.0D0
          do K=1,3
            qij=Q(K,I)-Q(K,J)
            if(qij > bl2) qij=qij-bl
            if(qij < -bl2) qij=qij+bl
            rsq=rsq+qij**2
          enddo
          if(rsq.GE.rc2) then
            QRC(I+(J-1)*N)=.FALSE.
          else
            QRC(I+(J-1)*N)=.TRUE.
          endif
          QRC(J+(I-1)*N)=QRC(I+(J-1)*N)
        enddo
      enddo
      return
end

