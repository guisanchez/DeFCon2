      subroutine write_plt(U,Z,celda,ncall,fin,ncp)

      implicit none

      integer nc,nn,nnx,nny,nvert,nvar,l,i,j,ncall
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar

      double precision Z(nc),U(nc,nvar),celda(nc,2)

      character*40 filename
      double precision xi,yi,dx,dy
      integer ncx, ncy,fin,ncp
      character*20 title
      integer node(4)

      title = 'titulo'

      ncx = nnx - 1
      ncy = nny - 1

      xi = 0.d0
      yi = 0.d0
      
      dx = abs(celda(1,1) - celda(2,1))
      dy = abs(celda(1,2) - celda(2,2))

C      if (fin.eq.0) then

C         if (ncall.eq.0) then
C            open(21,file='plototal.plt')  
C            write(21,*) 'TITLE=""'
C         else
C            open(21,file='plototal.plt',status='OLD',position='APPEND')
C         endif
C      else
C      if (fin.eq.1) then
C         open(21,file='plotfinal.vtk')
C         write(21,*) 'TITLE=""'
C      endif
      write(filename,'(a5,i4.4,a4)') 'plot_',ncall,'.vtk'
      write(*,*) ncall,filename
      open(19,file=filename)

      write(19,'(a26)')  '# vtk DataFile Version 2.0'
      write(19,'(a29)')   adjustl(title)
      write(19,'(a6)')  'ASCII '
      write(19,'(a26)')  'DATASET UNSTRUCTURED_GRID '

      write(19,'(a6,1X,i10.10,1X,a5)') 'POINTS', nnx*nny, 'float'

      do j = 1, ncy
         do i = 1, ncx
            l = ( j - 1 ) * ncx + i
            write(19,'(3(F10.3,1X))') celda(l,1) - dx*0.5, 
     +           celda(l,2) - dy*0.5 , 0.0
         enddo
         write(19,'(3(F10.3,1X))') celda(l,1) + dx*0.5, 
     +        celda(l,2) - dy*0.5 , 0.0
      enddo
      j = ncy
      do i = 1, ncx
         l = ( j - 1 ) * ncx + i
         write(19,'(3(F10.3,1X))') celda(l,1) - dx*0.5, 
     +        celda(l,2) + dy*0.5 , 0.0
      enddo
      write(19,'(3(F10.3,1X))') celda(l,1) + dx*0.5, celda(l,2) + dy*0.5 
     +     , 0.0
      
      write(19,*) 'CELLS', ncp, 5*ncp ! ncx*ncy, 5*ncx*ncy
      
      do j = 1, ncy
         do i = 1, ncx
            l = ( j - 1 ) * ncx + i
            if (Z(l).ge.-9990.d0) then
               node(1) = (j-1)*nnx + (i - 1)
               node(2) = (j-1)*nnx + i 
               node(3) = j*nnx + i 
               node(4) = j*nnx + i - 1
               write(19,*) 4, node
            endif
         enddo
      enddo
      
      write(19,*) 'CELL_TYPES ', ncp ! ncx*ncy
      
C      do j = 1, ncy
C         do i = 1, ncx
C            write(19,*) 5
C         enddo
C      enddo

      do j = 1, ncp
         write(19,*) 5
      enddo

      write(19,*) 'CELL_DATA', ncp ! ncx*ncy
      write(19,*) 'SCALARS D float'
      write(19,*) 'LOOKUP_TABLE default'

      do j = 1, ncy
         do i = 1, ncx
            l = ( j - 1 ) * ncx + i
            if (Z(l).gt.-9990.d0) write(19,*) Z(l) + U(l,1)
         enddo
      enddo

      write(19,*) 'SCALARS h float'
      write(19,*) 'LOOKUP_TABLE default'

      do j = 1, ncy
         do i = 1, ncx
            l = ( j - 1 ) * ncx + i
            if (Z(l).gt.-9990.d0) write(19,*) U(l,1)
         enddo
      enddo

      write(19,*) 'SCALARS u float'
      write(19,*) 'LOOKUP_TABLE default'

      do j = 1, ncy
         do i = 1, ncx
            l = ( j - 1 ) * ncx + i
            if (Z(l).gt.-9990.d0) write(19,*) U(l,2)
         enddo
      enddo

      write(19,*) 'SCALARS v float'
      write(19,*) 'LOOKUP_TABLE default'

      do j = 1, ncy
         do i = 1, ncx
            l = ( j - 1 ) * ncx + i
            if (Z(l).gt.-9990.d0) write(19,*) U(l,3)
         enddo
      enddo

      



C            if (i.lt.ncx) then
C               write(19,'(7(1X,D18.6))',ADVANCE='NO') celda(l,1), 
C     +              celda(l,2), Z(l), Z(l) + U(l,1), U(l,1:3)
C               write(21,'(7(1X,D18.6))',ADVANCE='NO') celda(l,1), 
C     +              celda(l,2), Z(l), Z(l) + U(l,1), U(l,1:3)
C            else
C               write(19,'(7(1X,D18.6))',ADVANCE='YES') celda(l,1),
C     +              celda(l,2), Z(l), Z(l) + U(l,1), U(l,1:3)
C               write(21,'(7(1X,D18.6))',ADVANCE='YES') celda(l,1), 
C     +              celda(l,2), Z(l), Z(l) + U(l,1), U(l,1:3)
C            endif
C         enddo
C      enddo

      close(19)
C      close(21)


      return 
      end
