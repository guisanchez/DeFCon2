      subroutine comparador(U,celda,sum2,malla)

      implicit none

      integer nc,nn,nnx,nny,nvert,nvar
      integer nvw,nhw
      double precision g,k,cfl
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar
      common /n_wall/ nvw,nhw
      common /genvars/ g,k,cfl
      double precision fric1, fric2,h1,h2
      common /friction/ fric1, fric2,h1,h2
      double precision celda(nc,2)
      double precision U(nc,nvar)
      double precision factor
C      parameter(factor=1.152d0)
      double precision ave, sum1, sum2, na1
      integer ncx2,ncy2
      double precision xmax2,xmin1,xmin2
      double precision ymax2,ymin1,ymin2
      double precision, dimension (:,:), allocatable :: fin
      double precision ar(500000)
      character*40 aux,malla,proc_file
      double precision dist, distmin,wet,x,y
      integer n
      double precision h, dif, diftot,dx
      integer i,j,l
      logical file_exists
      double precision, dimension (:), allocatable :: depth
      
      proc_file = 'h_obs_proc.dat'

      xmax2 = -1.d+6
      ymax2 = -1.d+6
      xmin2 = 1.d+6
      ymin2 = 1.d+6

C      open(11,file='err_proc')
C      read(11,*) i
C      close(11)
 
      inquire(file=proc_file,EXIST=file_exists)

      factor = 1.d0
C      if (i.eq.1) then
      if(file_exists) then
         open(2,file=proc_file,form='UNFORMATTED')
         read(2) sum1,ave
         factor = sum1/sum2
         factor = 1.d0
         open(unit=1,file=malla)
         read(1,*) aux
         read(1,*) aux, ncx2
         read(1,*) aux, ncy2
         read(1,*) aux, xmin1
         read(1,*) aux, ymin1
         read(1,*) aux, dx
         close(1)
         allocate(depth(ncy2*ncx2))
         diftot = 0.d0  
         sum2 = 0.d0
         na1 = 0.d0
         do i = 1, ncy2*ncx2
            read(2) dif
            depth(i) = dif
            h = U(i,1)
            sum2 = sum2 + h
            if (dif.lt.-9990.d0) dif = 0.d0
            if (h.lt.-9990.d0) h = 0.d0
            if (dif.gt.0.d0) na1 = na1 + (dif - ave)*(dif - ave)
            dif = dif - h*factor
            
            diftot = diftot + dif * dif
            
         enddo
         
         diftot = sqrt(diftot)
         
         if (na1.gt.0.d0) then
            na1 = 1.d0 - diftot*diftot/na1
            open(111,file='Nash_Sutcliffe',POSITION='APPEND')
            write(111,*) 1.d0 + 1.d-10 - na1,fric1,fric2
            close(111)
            write(*,*) 1.d0 + 1.d-10 - na1
            
         endif
      
C         open(1,file='plotcomp.plt')
C         write(1,*) 'TITLE=""'
C         write(1,*) 'VARIABLES="X","Y","H1","H2"'
C         write(1,*) 'ZONE I=',ncx2,', J=',ncy2,
C     +        ', DATAPACKING=POINT'
C         do j = 1, ncy2
C            do i = 1, ncx2
C               l = ( j - 1 ) * ncx2 + i
C               write(1,'(7(1X,D18.6))',ADVANCE='YES') celda(l,1), 
C     +              celda(l,2),U(l,1),depth(l)
C            enddo
C         enddo
C         
C         close(1)

      else

         open(2,file='h_obsfinalrot.tec')
         read(2,*) aux
         read(2,*) aux
         read(2,*) aux
         ncx2 = 267
         ncy2 = 260
         n = ncx2*ncy2
         allocate(fin(n,3))
         j = 0
         wet = 0.d0
         sum1 = 0.d0
         do i = 1, ncy2
            read(2,*) ar(1:ncx2*3)
            do j = 1, ncx2
               l = (i-1)*ncx2 + j
               fin(l,1) = ar((j-1) * 3 + 1)
               fin(l,2) = ar((j-1) * 3 + 2)
               fin(l,3) = ar(j*3)
               if (fin(l,3).gt.0.d0) then
                  sum1 = sum1 + fin(l,3)
                  wet = wet + 1.d0
               endif
               if (fin(l,1).gt.xmax2) xmax2 = fin(l,1)
               if (fin(l,1).lt.xmin2) xmin2 = fin(l,1)
               if (fin(l,2).gt.ymax2) ymax2 = fin(l,2)
               if (fin(l,2).lt.ymin2) ymin2 = fin(l,2)
            enddo
         enddo
         if (wet.gt.0.d0) then
            ave = sum1 / wet
         endif
!     !      write(*,*) sum1

         open(unit=1,file=malla)
         read(1,*) aux
         read(1,*) aux, ncx2
         read(1,*) aux, ncy2
         read(1,*) aux, xmin1
         read(1,*) aux, ymin1
         read(1,*) aux, dx
         close(1)
         
         diftot = 0.d0  
         sum2 = 0.d0
         do i = 1, ncy2
            do j = 1, ncx2
               dif = 0.d0
               x = xmin1 + real(j-1)*dx
               y = ymin1 + real(i-1)*dx
               h = U((i-1)*ncx2+j,1)
               sum2 = sum2 + h
               if (((x-xmax2)*(x-xmin2).lt.0.d0).and.
     +              ((y-ymax2)*(y-ymin2).lt.0.d0)) then
                  distmin = 1.d+6
                  do  l = 1, n
                     dist = (x - fin(l,1))*(x - fin(l,1))
                     dist = dist + (y - fin(l,2))*(y - fin(l,2))
                     dist = sqrt(dist)
                     if (dist.lt.distmin) then
                        distmin = dist
                        dif = fin(l,3)
                     endif
                  enddo
                  if (dif.lt.-9990.d0) dif = 0.d0
                  if (h.lt.-9990.d0) h = 0.d0
C     if (distmin.gt.1.d0) write(*,*) distmin
               endif
               dif = dif - h*factor
C     write(10,*) x, y, dif
               diftot = diftot + dif * dif
            enddo
         enddo
!     !      write(*,*) sum2
         diftot = sqrt(diftot)
         na1 = 0.d0
         do l = 1, n
            if (fin(l,3).gt.0.d0) then
               na1 = na1 + (fin(l,3) - ave)*(fin(l,3) - ave)
            endif
         enddo
         if (na1.gt.0.d0) then
            na1 = 1.d0 - diftot*diftot/na1
            open(111,file='Nash_Sutcliffe',POSITION='APPEND')
            write(111,*) 1.d0 + 1.d-10 - na1,fric1,fric2
            close(111)
            write(*,*) 1.d0 + 1.d-10 - na1
!     !         write(*,*) 'E=',na1
         endif
      endif
      return

      end
