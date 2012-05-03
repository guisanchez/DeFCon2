
      subroutine read_mesh2D(node,celda,hw,vw,hwn,vwn,cellw,cellnod,Z,U
     +     ,ncp)

      implicit none

      integer nn                ! Number of nodes
      integer nnx,nny           ! Number of nodes in x and y directions, resp.
      integer nc                ! Number of cells = (nnx-1)*(nny-1)
      integer nvert             ! Number of vert. (4)
      integer nvar              ! Number of cons. variables
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar
      integer nhw,nvw           ! num. horizontal walls, num. vertical walls
      common /n_wall/ nvw, nhw
      integer nw
      integer i,j,l,au1
      integer ncx, ncy,k

      double precision node(nn,2)
      double precision celda(nc,2) ! Centro de masas de la celda
      integer cellnod(nc,nvert)
      integer cellw(nc,nvert) ! Las paredes de la celda
C     cellw(i,1:2) = paredes horizontales; cellw(i,3:4) = paredes verticales
C      double precision cellvar(nc,nvar)
      double precision Z(nc),U(nc,nvar)

      integer hw(nhw,2),vw(nvw,2),hwn(nhw,2),vwn(nvw,2)
      integer pared(nhw+nvw,2)
      character *40 achar
C     hw = horizontal wall, vw = vertical wall 
C     hw(i,1:2) are cells over and below horizontal wall 'i' 
C     hwn = horizontal wall nodes, vwn = vertical wall nodes
C     vwn(j,1:2) are nodes that determine vertical wall 'j'
      double precision x_min,x_max,y_min,y_max,dx,dy,eps,xn0,yn0,ynl,xnl
      double precision pi
      parameter (pi = 3.141592653)
      integer nodata,preproc
      integer piz,ps,pld,pli
      integer ncp
      
      open(11,file='preproc')
      
      read(11,*) preproc
      
      close(11)
      if (preproc.eq.1) then
C         write(*,*) 'Reading proccessed mesh from binary file'
         open(10,file='proc_mesh.dat',form='UNFORMATTED')
C         write(*,*) 'reading node positions'
         do i = 1, nn
            read(10) node(i,1:2)
         enddo
C         write(*,*) 'reading cell centers'
         do i = 1, nc
            read(10) celda(i,1:2)
         enddo
C         write(*,*) 'reading horizontal wall cells'
         do i = 1, nhw
            read(10) hw(i,1:2)
         enddo
C         write(*,*) 'reading vertical wall cells'
         do i = 1, nvw
            read(10) vw(i,1:2)
         enddo
C         write(*,*) 'reading horizontal wall nodes'
         do i = 1, nhw
            read(10) hwn(i,1:2)
         enddo
C         write(*,*) 'reading vertical wall nodes'         
         do i = 1, nvw
            read(10) vwn(i,1:2)
         enddo
C         write(*,*) 'reading cell walls' 
         do i = 1, nc
            read(10) cellw(i,1:nvert)
         enddo
C         write(*,*) 'reading cell nodes' 
         do i = 1, nc
            read(10) cellnod(i,1:nvert)
         enddo
C         write(*,*) 'reading cell Z' 
         do i = 1, nc
            read(10) Z(i)
         enddo
         close(10)
C         write(*,*) 'opening init0.txt'
         open(2,file='init0.txt')
         do i = 1, 7
            read(2,*)
         enddo
         ncx = nnx - 1 
         ncy = nny - 1
         do j = 1, ncy
C            write(*,*) j,ncy
            k = (j-1)*ncx + 1
            read(2,*) U(k:k+ncx-1,1)
         enddo
         close(2)
      else
C         write(*,*) 'Processing mesh...'
         nw = nhw + nvw 

         x_min = 1.d+6
         y_min = 1.d+6
         x_max = -1.d+6
         y_max = -1.d+6

C     Nodos y paredes en celda rectangular C_i
C
C             W_i2
C               |
C               v
C         N_i4-----N_i3
C          |        |                    ^
C  W_i3->  |   C_i  |  <-W_i4            |
C          |        |                  y |
C         N_i1-----N_i2                  |------>
C              ^                            x
C              |
C            W_i1
C
C     Orden nodos: N_i1 < N_i2 < N_i3 < N_i4
C     Paredes: Separación Paredes || x ,  ||y
C     Horizontales: W_i1 < W_i2
C     Verticales: W_i3 < W_i4
C
C      write(*,*) nc,nn,nnx,nny,nvert,nvar
C      write(*,*) nvw, nhw
C      stop

         open(unit=1,file='malla0.txt')
         open(unit=2,file='init0.txt')
C     Read number of vert
         read(1,*) achar,nvert
         read(2,*)
C     read number of y-cells and x-cells
         read(1,*) achar, ncx
         read(1,*) achar, ncy
         read(2,*)
         read(2,*)
C     number of y-nodes and x-nodes
         nnx = ncx + 1
         nny = ncy + 1
C     read (x0,y0) origin
         read(1,*) achar, x_min
         read(1,*) achar, y_min
         read(2,*)
         read(2,*)
C     read dx, dy
         read(1,*) achar, dx
C     write(*,*) dx
         read(2,*)
         dy = dx
C     Read nodata value
         read(1,*) achar,nodata
         read(2,*)
C     Number of cells
         nc = (nnx-1)*(nny -1)
C     Number of nodes
         nn = nnx * nny
C     xmax, ymax (cell centers)
         x_max = x_min + dx * (ncx-1)
         y_max = y_min + dy * (ncy-1) 
C     xn0, yn0:: min x,y for nodes
         xn0 = x_min - dx*0.5d0
         yn0 = y_min - dy*0.5d0
C     xnl, ynl :: max x,y for nodes
         xnl = x_max + dx*0.5d0
         ynl = y_max + dy*0.5d0
         
C     Node positions
C     write(*,*) 'Node positions'
         do j = 1, nny
            do i = 1, nnx
               k = (j-1)*nnx + i
               node(k,1) = xn0 + real(i-1) * dx
               node(k,2) = yn0 + real(j-1) * dy
            enddo
         enddo
         
         
C      k = 1
C      write(*,*) 'node',k,':'
C      write(*,*) node(k,1:2)
c      k = 114577
c      write(*,*) 'node',k,':'
c      write(*,*) node(k,1:2)
c      k = 114578
c      write(*,*) 'node',k,':'
c      write(*,*) node(k,1:2)
c      k = 114918
c      write(*,*) 'node',k,':'
c      write(*,*) node(k,1:2)
C      stop

C     Cells
         j = 0
C      write(*,*) 'Cell nodes'
         do i = 1, nn - nnx       
            if (mod(i,nnx).ne.0) then
               j = j + 1
C     write(*,*) j
               cellnod(j,1) = i
               cellnod(j,2) = i + 1
               cellnod(j,3) = i + 1 + nnx
               cellnod(j,4) = i + nnx
C            if (j.eq.1) then
C               write(*,*) '__________'
C               write(*,*) 'cell',j,':'
C               write(*,*) nc,j,cellnod(1,1:4)
C            endif
               celda(j,1) = 0.25*(node(cellnod(j,1),1) + 
     +              node(cellnod(j,2),1) + node(cellnod(j,3),1) + 
     +              node(cellnod(j,4),1))
               celda(j,2) = 0.25*(node(cellnod(j,1),2) + 
     +              node(cellnod(j,2),2) + node(cellnod(j,3),2) + 
     +              node(cellnod(j,4),2))
            endif
         enddo
      
C      j = 1
C      write(*,*) 'cell',j,':'
C      write(*,*) 'nodes',cellnod(j,1:4)
C      stop
      
C
C     READ NODE DATA
C
C      write(*,*) 'Nodes data'
         do j = 1, ncy
            k = (j-1)*ncx + 1
            read(1,*) Z(k:k+ncx-1)
            read(2,*) U(k:k+ncx-1,1)
C         write(*,*) j,k,':',k+ncx,';',nc
C         do i = k, k + ncx
C            write(*,*) i,Z(i)
C         enddo
         enddo
         close(1)
         close(2)
      
         nvw = nnx*(nny-1)
         nhw = (nnx-1)*nny
      

C                                                                    C
C     NODOS QUE DEFINEN UNA PARED Y CELDAS QUE DEFINEN UNA PARED     C
C                                                                    C

C     Paredes horizontales
C      write(*,*) 'horizontal walls'
         do i = 1, nc
C     Por decreto: la pared horizontal i está bajo la celda i
C     i.e., sus nodos son los nodos inferiores de la celda i
            hwn(i,1:2) = cellnod(i,1:2)  
C     además, la celda sobre la pared i es la celda i
            hw(i,2) = i
         enddo
      
         i = nc
         eps = dy/1000.d0

         au1 = 0

         do j = 1, nc
C     Celdas de la fila superior
            if (abs(celda(j,2)-y_max).lt.eps) then
               i = i + 1
C     Estamos en la fila superior: los nodos de la pared i
C     son los nodos superiores de la celda j
               hwn(i,1) = cellnod(j,4)
               hwn(i,2) = cellnod(j,3)
C     además la celda bajo la pared i es la celda j
               hw(i,1) = j
C     y la celda sobre la pared i no existe
               hw(i,2) = 0
            endif
            if (i.eq.nhw) exit
         enddo

C     El otro vecino de cada pared
         do i = 1, nhw
C     Paredes de la fila inferior
            if (abs(node(hwn(i,1),2)-yn0).lt.eps) then
               hw(i,1) = 0      ! estas no tienen vecino debajo
C     A continuación las demás paredes que no son fila superior
            else if (abs(node(hwn(i,1),2)-ynl).gt.eps) then
C     Veamos cual es la celda que tiene debajo la pared i
               do j = 1, nc
                  if (((cellnod(j,3).eq.hwn(i,2)).and.
     +                 (cellnod(j,4).eq.hwn(i,1))).or.
     +                 ((cellnod(j,3).eq.hwn(i,1)).and.
     +                 (cellnod(j,4).eq.hwn(i,2)))) then
                     hw(i,1) = j
                     exit
                  endif
               enddo
            endif
            
         
         enddo

C     Paredes verticales
C      write(*,*) 'vertical walls'
         do i = 1, nc
            l = i + nhw
C     Por decreto: la pared vertical i está a la izquierda de la celda i
C     i.e., sus nodos son los nodos izquierdos de la celda i
            vwn(i,1) = cellnod(i,1)
            vwn(i,2) = cellnod(i,4)
C     por tanto, la celda a la derecha la pared i es la celda i
            vw(i,2) = i

         enddo
         eps = dx/1000.d0         
         i = nc
         do j = 1, nc
C     Celdas del lateral derecho
            if (abs(celda(j,1)-x_max).lt.eps) then
               i = i + 1
C     Estamos en el lateral derecho: los nodos de la pared i
C     son los nodos derechos de la celda j
               vwn(i,1) = cellnod(j,2)
               vwn(i,2) = cellnod(j,3)
C     además la celda a la izquierda de la pared i es la celda j
               vw(i,1) = j
C     y la celda a la derecha de la pared i no existe
               vw(i,2) = 0
            endif
            if (i.eq.nvw) exit
         enddo

C     El otro vecino de cada pared
         do i = 1, nvw
C     Paredes del lateral izquierdo 
            if (abs(node(vwn(i,1),1)-xn0).lt.eps) then
C            if (i.eq.80) write(*,*) i,'lat. izdo'
               vw(i,1) = 0      ! estas paredes no tienen vecino izquierdo
C     A continuación las demás paredes que no son lateral del dominio
C         else if (abs(node(cellnod(i,2),1)-x_max).gt.eps) then
            else if (abs(node(vwn(i,2),1)-xnl).gt.eps) then
C            if (i.eq.80) write(*,*) i,'no lateral'
C     Veamos cual es la celda que tiene a la izquierda la pared i
               do j = 1, nc
                  if (((cellnod(j,2).eq.vwn(i,2)).and.
     +                 (cellnod(j,3).eq.vwn(i,1))).or.
     +                 ((cellnod(j,2).eq.vwn(i,1)).and.
     +                 (cellnod(j,3).eq.vwn(i,2)))) then
                     vw(i,1) = j
                     exit
                  endif
               enddo
            endif
         enddo
         


         pared(1:nhw,1:2) = hw(1:nhw,1:2)
         pared(1+nhw:nw,1:2) = vw(1:nvw,1:2)

C     Sabemos las celdas que definen una pared
C     A continuación: las paredes que definen una celda
C      write(*,*) 'cell walls'
         cellw = 0
         do j = 1, nc        
            do i = 1, nhw
C     celda j. 
C     Si es la que la pared i tiene debajo, entonces i es la 
C     pared horizontal superior de la celda j
C     Si es la que la pared i tiene encima, entonces i es la
C     pared horizontal inferior de la celda j
               if (hw(i,1).eq.j) cellw(j,2) = i
               if (hw(i,2).eq.j) cellw(j,1) = i
               if (cellw(j,1)*cellw(j,2).ne.0) exit
            enddo
            do i = 1, nvw
C     celda j.
C     Si es la que la pared i tiene a la izquierda, entonce i es la
C     pared lateral derecha de la celda j
C     Si es la que la pared i tiene a su derecha, entonces i es la 
C     pared lateral izquierda de la celda j
               if (vw(i,1).eq.j) cellw(j,3) = i
               if (vw(i,2).eq.j) cellw(j,4) = i
               if (cellw(j,3)*cellw(j,4).ne.0) exit
            enddo
         enddo
      endif




C     Eliminación de celdas que no forman parte del dominio de interés
C     ELIMINACIÓN CELDAS CON U(i,1) <= -9999
C      write(*,*) 'WARNING' 
C      write(*,*) 'Cells with Z <= -9999 will be "discarded"'
C      write(*,*) 'Negative h cells: h -> -h'
C      write(*,*) 'flat bed'
c      write(*,*) 'Cells with h <= 0 will be set to 0'
      ncp = 0
      do i = 1, nc
         U(i,2:nvar) = 0.d0
         if (U(i,1).lt.0.d0) then
            if (U(i,1).gt.-9998.d0) then
               U(i,1) = -U(i,1)
            else
               U(i,1) = 0.d0
            endif
         endif
         ncp = ncp + 1
         if (Z(i).lt.-9998.d0) then
            ncp = ncp -1
C            Z(i) = 0.d0
C     Pared inferior
            piz = cellw(i,1)
C     Pared superior
            ps = cellw(i,2)
C     Pared lateral izquierda
            pli = cellw(i,4)
C     Pared lateral derecha
            pld = cellw(i,3)
            hw(piz,2) = 0
            hw(ps,1) = 0
            vw(pli,2) = 0
            vw(pld,1) = 0
         endif
      enddo
      return
      
      end
      




 
