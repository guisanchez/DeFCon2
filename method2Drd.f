
      subroutine method2D(node,hw,vw,hwn,vwn,cellw,U
     +     ,dUp,dUn,dt,nstep)

      implicit none

      integer nc,nn,nnx,nny,nvert,nvar,nstep
      integer nvw,nhw
      double precision g,k,cfl
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar
      common /n_wall/ nvw,nhw
      common /genvars/ g,k,cfl

      integer hw(nhw,2),vw(nvw,2),hwn(nhw,2),vwn(nvw,2)
      double precision dt
      
      integer i

      double precision node(nn,2)
      integer cellw(nc,nvert)
      double precision U(nc,nvar),Un(nc,nvar)
      double precision Un2(nc,nvar),dU(nvar)

      double precision dUp(nhw+nvw,nvar,2),dUn(nhw+nvw,nvar,2)

      double precision dx,dy,Ai
      integer pi,ps,pli,pld ! Pared inferior, superior,lateral izda, lat. dcha
      integer noflux
      double precision suma

      suma = 0.d0

C      do i = 1, nc
C         suma = suma + U(i,1)
C      enddo

      do i = 1, nc   ! loop over any cell

         dU(1:nvar) = 0.d0
C                            C
C     Área de la celda i     C
C                            C

C     Pared inferior (horizontal)
         pi = cellw(i,1)
C     Nodos que delimitan pared pi: node(hwn(pi,1),1:2) <->node(hwn(pi,2),1:2)
         dx = abs(node(hwn(pi,1),1) - node(hwn(pi,2),1))

C     Pared lateral izquierda (vertical)
         pli = cellw(i,4) !+ nhw
C     Nodos que delimitan pared pl: node(vwn(pli,1),1:2) <->node(vwn(pli,2),1:2)
         dy = abs(node(vwn(pli,1),2) - node(vwn(pli,2),2))

C     Área
         Ai = dx*dy

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CALCULANDO PREDICTOR    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                             C
C     Contribuciones paredes horizontales     C
C                                             C
         
C     Pared inferior
         pi = cellw(i,1)
         dU(1:nvar) = dU(1:nvar) + dUp(pi,1:nvar,1)
C     Pared superior
         ps = cellw(i,2)
         dU(1:nvar) = dU(1:nvar) + dUn(ps,1:nvar,1)
C     Pared lateral izquierda
         pli = cellw(i,4) + nhw
         dU(1:nvar) = dU(1:nvar) + dUp(pli,1:nvar,1)
C     Pared lateral derecha
         pld = cellw(i,3) + nhw
         dU(1:nvar) = dU(1:nvar) + dUn(pld,1:nvar,1)

         Un(i,1:nvar) = U(i,1:nvar) + dU(1:nvar)*dt/Ai
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CALCULANDO CORRECTOR     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         dU(1:nvar) = 0.d0
         
C     Pared inferior
         dU(1:nvar) = dU(1:nvar) + dUp(pi,1:nvar,2)
C     Pared superior
         dU(1:nvar) = dU(1:nvar) + dUn(ps,1:nvar,2)
C     Pared lateral izquierda
         dU(1:nvar) = dU(1:nvar) + dUp(pli,1:nvar,2)
C     Pared lateral derecha
         dU(1:nvar) = dU(1:nvar) + dUn(pld,1:nvar,2)
     
         Un2(i,1:nvar) = Un(i,1:nvar) + dU(1:nvar)*dt/Ai

         if (Un2(i,2)*Un(i,2).lt.0.d0) then
            Un2(i,2) = 0.d0
         endif
         if (Un2(i,3)*Un(i,3).lt.0.d0) then
            Un2(i,3) = 0.d0
         endif

         U(i,1:nvar) = Un2(i,1:nvar)
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CONDICIONES DE CONTORNO     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Celdas en el borde superior o inferior
C     noflux = id_celda_inferior*id_celda_superior
C     Si la celda i está en el borde, alguna de ellas tendrá id = 0
C     si noflux = 0 -> U(i,3) = 0 (Si la celda está arriba o abajo, qy=0)
         noflux = hw(cellw(i,1),1)*hw(cellw(i,2),2)
         if (noflux.eq.0) then
            U(i,3) = 0.d0
         endif
C     Celdas en el borde derecho o izquierdo
C     noflux = *id_celda_izda*id_celda_dcha         
C     Si la celda i está en el borde, alguna de ellas tendrá id = 0
C     si noflux = 0 -> U(i,3) = 0 (Si la celda está a dcha o izda, qx=0)
         noflux = vw(cellw(i,4),1)*vw(cellw(i,3),2)
         if (noflux.eq.0) then
            U(i,2) = 0.d0
         endif

      enddo

      suma = 0.d0

C      do i = 1, nc
C         suma = suma + U(i,1)
C      enddo

      call redis(hw,vw,cellw,U)
      
      suma = 0.d0

C      do i = 1, nc
C         suma = suma + U(i,1)
C      enddo
      
      return

      end

      subroutine redis(hw,vw,cellw,U)

C                                                        C
C     Por una redistribución de la riqueza más justa     C
C                                                        C

      implicit none

      integer nc,nn,nnx,nny,nvert,nvar
      integer nvw,nhw
      double precision g,k,cfl
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar
      common /n_wall/ nvw,nhw
      common /genvars/ g,k,cfl

      integer hw(nhw,2),vw(nvw,2)

      
      integer i


      integer cellw(nc,nvert)
      double precision U(nc,nvar)
      integer pi,ps,pli,pld ! Pared inferior, superior,lateral izda, lat. dcha
      integer ci,cs,cli,cld ! Celda inferior, superior,lateral izda, lat. dcha
      double precision hpos(4)
      double precision h_ne,theta

      do i = 1, nc
C     Pared inferior
         pi = cellw(i,1)
C     Pared superior
         ps = cellw(i,2)
C     Pared lateral izquierda
         pli = cellw(i,4) + nhw
C     Pared lateral derecha
         pld = cellw(i,3) + nhw
C     Celda inferior
         ci = hw(cellw(i,1),1)
C     Celda superior
         cs = hw(cellw(i,2),2)
C     Celda lateral izquierda
         cli = vw(cellw(i,4),1)
C     Celda lateral derecha
         cld = vw(cellw(i,3),2)
         
         hpos (1:4) = 0.d0
         h_ne = 0.d0
         if (U(i,1).lt.0.d0) then  ! ¡¡Calado negativo !!
            if ((ci.ne.0).and.(U(ci,1).gt.0.d0)) then
               hpos(1) = 1.d0
               h_ne = h_ne + U(ci,1)
            endif
            if ((cs.ne.0).and.(U(cs,1).gt.0.d0)) then
               hpos(2) = 1.d0
               h_ne = h_ne + U(cs,1)
            endif
            if ((cli.ne.0).and.(U(cli,1).gt.0.d0)) then
               hpos(3) = 1.d0
               h_ne = h_ne + U(cli,1)
            endif
            if ((cld.ne.0).and.(U(cld,1).gt.0.d0)) then
               hpos(4) = 1.d0
               h_ne = h_ne + U(cld,1)
            endif

            if (h_ne.lt.-U(i,1)) then
               goto 1000
            else
C     Reparto equitativo
               theta = (h_ne + U(i,1))/h_ne
               if ((cld.ne.0).and.(U(cld,1).gt.0.d0)) 
     +              U(cld,1) = U(cld,1) * theta
               if ((cli.ne.0).and.(U(cli,1).gt.0.d0)) 
     +              U(cli,1) = U(cli,1) * theta
               if ((cs.ne.0).and.(U(cs,1).gt.0.d0)) 
     +              U(cs,1) = U(cs,1) * theta
               if ((ci.ne.0).and.(U(ci,1).gt.0.d0)) 
     +              U(ci,1) = U(ci,1) * theta
               U(i,1) = 0.d0
            endif
         endif
 1000    continue         
      enddo

      return

      end
