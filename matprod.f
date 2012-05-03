      subroutine matprod(m1,m2,m3,MP,n)
      
      implicit none
      
      integer n,i,j,k
C      parameter (n=2)
      double precision m1(n,n), m2(n,n), m3(n,n),MP(n,n),m_aux(n,n)

      m_aux = 0.0
      do i = 1, n
         do j = 1, n
            do k = 1, n
               m_aux(i,j) = m_aux(i,j) + m1(i,k)*m2(k,j)
            enddo
         enddo
      enddo

      MP = 0.0
      do i = 1, n
         do j = 1, n
            do k = 1, n
               MP(i,j) = MP(i,j) + m_aux(i,k)*m3(k,j)
            enddo
         enddo
      enddo
C      write(96,*) m1(1,1:2),m2(1,1:2),m3(1,1:2),MP(1,1:2)
C      write(96,*) m1(2,1:2),m2(2,1:2),m3(2,1:2),MP(2,1:2)
      return
      end
