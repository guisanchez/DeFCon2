      program preproc_comparador

      implicit none

      integer nc,nn,nnx,nny,nvert,nvar
      integer nvw,nhw
      double precision g,k,cfl
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar
      common /n_wall/ nvw,nhw
      common /genvars/ g,k,cfl
      double precision fric1, fric2,h1,h2
      common /friction/ fric1, fric2,h1,h2
      double precision factor
      parameter(factor=1.152d0)
      double precision ave, sum1, sum2, wet
      integer n, ncx2,ncy2
      double precision xmax2,xmin1,xmin2
      double precision ymax2,ymin1,ymin2
      double precision, dimension (:,:), allocatable :: fin
      character*40 aux
      double precision dist, distmin, x, y, dif, diftot,dx
      integer i,j,l
      double precision ar(500000)
      

      xmax2 = -1.d+6
      ymax2 = -1.d+6
      xmin2 = 1.d+6
      ymin2 = 1.d+6

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

      open(unit=1,file='malla0.txt')
      read(1,*) aux
      read(1,*) aux, ncx2
      read(1,*) aux, ncy2
      read(1,*) aux, xmin1
      read(1,*) aux, ymin1
      read(1,*) aux, dx
      close(1)
      open(2,file='h_obs_proc.dat',form='UNFORMATTED')
      write(2) sum1,ave
      diftot = 0.d0  
      sum2 = 0.d0
      do i = 1, ncy2
         do j = 1, ncx2
            dif = 0.d0
            x = xmin1 + real(j-1)*dx
            y = ymin1 + real(i-1)*dx
            if (((x-xmax2)*(x-xmin2).lt.0.d0).and.
     +           ((y-ymax2)*(y-ymin2).lt.0.d0)) then
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

               if (distmin.lt.dx) then
                  write(2) dif
               else
                  write(2) 0.d0
               endif
            else
               write(2) 0.d0
            endif
         enddo
      enddo
      close(2)

      open(10,file='err_proc')
      write(10,*) 1
      close(10)

      end
