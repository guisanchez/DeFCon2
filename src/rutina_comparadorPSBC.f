      subroutine comparador(U,celda,sum2,malla,h_obs,outns)

      implicit none

      integer nc,nn,nnx,nny,nvert,nvar
      integer nvw,nhw
      double precision g,k,cfl
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar
      common /n_wall/ nvw,nhw
      common /genvars/ g,k,cfl
      double precision fric1, fric2,fric3,fr_hmin,fr_qmin,rho,h1,h2
      integer frstyle
      common /friction/ fric1, fric2, fric3,frstyle,fr_hmin,fr_qmin,rho

      double precision celda(nc,2)
      double precision U(nc,nvar)
      double precision factor
C      parameter(factor=1.152d0)
      double precision ave, sum1, sum2, na1
      integer ncx2,ncy2,ncx1,ncy1
      double precision xmax2,xmin1,xmin2
      double precision ymax2,ymin1,ymin2
      double precision, dimension (:,:), allocatable :: fin
      double precision ar(500000)
      character*40 aux,malla,proc_file
      double precision dist, distmin,wet,x,y
      integer n
      double precision h, dif, diftot,dx
      integer i,j,l,l1l
      logical file_exists
      double precision, dimension (:), allocatable :: depth
      character*40 h_obs
      integer p_err, p_mesh
      common /preproc/ p_err,p_mesh
      character*40 outns

      proc_file = 'h_obs_proc.dat'

      xmax2 = -1.d+6
      ymax2 = -1.d+6
      xmin2 = 1.d+6
      ymin2 = 1.d+6

 
      inquire(file=proc_file,EXIST=file_exists)

      factor = 1.d0

      if(file_exists) then
C         write(*,*) 'Error preprocesado'
C         call system('date')
         open(2,file=proc_file,form='UNFORMATTED')
         read(2) sum1,ave
         factor = sum1/sum2
C         write(*,*) 'factor = ',factor
C         factor = 1.d0
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
C            open(111,file='Nash_Sutcliffe',POSITION='APPEND')
			open(111,file=outns)
            if (frstyle.eq.1) then
               write(111,*) 1.d0 + 1.d-10 - na1,fric1,fric2,k,'V'
               close(111)
            else if (frstyle.eq.2) then
               write(111,*) 1.d0 + 1.d-10 - na1,rho,fric1,fric2,k,'B' 
               close(111)
            else
               write(111,*) 1.d0+1.d-10-na1,rho,fric1,fric2,fric3,k,'C'
               close(111)
            endif
            write(*,*) 1.d0 + 1.d-10 - na1
            
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      else
         
C         write(*,*) 'Procesando el error'
C         call system('date')
C         write(*,*) h_obs
         if (p_err.eq.1) then
            open(22,file = 'h_obs_proc.dat',form='UNFORMATTED')
         endif
         open(2,file=h_obs)
         read(2,*) aux

         read(2,*) aux,ncx2

         read(2,*) aux,ncy2
         read(2,*) aux,xmin2

         read(2,*) aux,ymin2

         read(2,*) aux,dx

         read(2,*) aux
         n = ncx2*ncy2
         allocate(fin(n,3))
         xmax2 = xmin2 + real(ncx2-1)*dx
         ymax2 = ymin2 + real(ncy2-1)*dx
         j = 0
         wet = 0.d0
         sum1 = 0.d0
         do i = 1, ncy2

            read(2,*) ar(1:ncx2)
            do j=1, ncx2
               if (ar(j).gt.0.d0) then
                  wet = wet + 1.d0
                  sum1 = sum1 + ar(j)
               endif
               l = (i-1)*ncx2 + j
               fin(l,3) = ar(j)  ! la h
               fin(l,1) = xmin2 + real(j-1)*dx  ! la x
               fin(l,2) = ymin2 + real(i-1)*dx  ! la y
            enddo

         enddo
         factor = sum1/sum2
C         write(*,*) 'factor = ',factor
         close(2)
         if (wet.gt.0.d0) then
            ave = sum1 / wet
         endif

         if (p_err.eq.1) then
            write(22) sum1, ave
         endif

         open(unit=1,file=malla)
         read(1,*) aux
         read(1,*) aux, ncx1
         read(1,*) aux, ncy1
         read(1,*) aux, xmin1
         read(1,*) aux, ymin1
         read(1,*) aux, dx
         close(1)
         
         diftot = 0.d0  
         sum2 = 0.d0
         do i = 1, ncy1
            do j = 1, ncx1
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
                        l1l = l
                     endif
                  enddo

                  if (dif.lt.-9990.d0) dif = 0.d0
                  if (h.lt.-9990.d0) h = 0.d0
                  
                  if (distmin.lt.dx) then
                     if (p_err.eq.1)   write(22) dif
                  else
                     if (p_err.eq.1) write(22) 0.d0
                  endif
               else
                  if (p_err.eq.1)  write(22) 0.d0
               endif
               dif = dif - h*factor
               diftot = diftot + dif * dif
            enddo
         enddo

         diftot = sqrt(diftot)
         na1 = 0.d0
         do l = 1, n
            if (fin(l,3).gt.0.d0) then
               na1 = na1 + (fin(l,3) - ave)*(fin(l,3) - ave)
            endif
         enddo
         if (na1.gt.0.d0) then
            na1 = 1.d0 - diftot*diftot/na1
C            open(111,file='Nash_Sutcliffe',POSITION='APPEND')
			open(111,file=outns)
            write(111,*) 1.d0 + 1.d-10 - na1,fric1,fric2,k
            close(111)
            write(*,*) 1.d0 + 1.d-10 - na1

         endif
      endif
C      call system('date')
      return

      end
