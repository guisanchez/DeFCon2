C     VERSION DE OPTIMIZACION DE PARAMETROS DE ENTRADA
C     ELIMINAR LOS C!C SI SE QUIERE ESCRIBIR FICHEROS DE SALIDA
C     ELIMINAR LOS !!C SI SE QUIERE VER LA SALIDA ESTÁNDAR

!> @file main.f
!< @param nc
!! @todo probando doxygen

      program Debris_2D

      implicit none

      integer i,j,l,lc,rc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Number of dimensions and number of cell variables
      integer ndim, nvar 
C      parameter (ndim = 2, nvar = 3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Variables
C     Nodes and cells
      integer nn                ! Number of nodes
      integer nnx,nny           ! Number of nodes in x and y directions, resp.
      integer nc                ! Number of cells = (nnx-1)*(nny-1)
      integer nvert             ! Number of vert. (4)
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar
      integer ncx, ncy

      integer nvw, nhw, nw      ! Number of vertical, horizontal and total walls
      common /n_wall/ nvw, nhw
      integer, dimension (:,:),allocatable ::  hw,vw,hwn,vwn
C     hw = horizontal wall, vw = vertical wall 
C     hw(i,1:2) are cells over and below horizontal wall 'i' 
C     hwn = horizontal wall nodes, vwn = vertical wall nodes
C     vwn(j,1:2) are nodes that determine vertical wall 'j'
      integer, dimension (:,:), allocatable :: cellw, cellnod
C     cellw(1:nc,1:4) = Las cuatro paredes de cada celda
C     cellw(i,1:2) = paredes horizontales; cellw(i,3:4) = paredes verticales
C     cellnod(1:nc,1:nvert) Los cuatro nodos de cada celda

      double precision, dimension (:,:), allocatable :: node, celda
      
C     node(k,1) = x position of k node, node(k,2) = y pos.
      double precision, dimension (:), allocatable :: Z ! Z(i) = bed height of i cell
      double precision, dimension (:,:), allocatable :: angle 
C     Slope: angle(k,1) = angle in x direction. angle(k,2) = angle in y direct.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Cell variables
      double precision, dimension (:,:), allocatable :: U , F1, G1
C     U -> variables de cada celda.
C     U(i,1) = h de celda i.
C     U(i,2) = h*u de la celda i.
C     U(i,3) = h*v de la celda i.
C     F1 -> flujo en x de cada celda.
C     F1(i,1) = h*u de celda i.
C     F1(i,2) = h*u*u de celda i.
C     F1(i,3) = h*u*v de celda i.
C     G1 -> flujo en y de cada celda.
C     G1(i,1) = h*v de celda i.
C     G1(i,2) = h*v*u de celda i.
C     G1(i,3) = h*v*v de celda i.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Variables in matrix form
      double precision, dimension (:,:,:), allocatable :: dUp,dUn
C     x = p,n
C     dUx(i,j,1) -> dUx en nodo i, componente j, predicho sin rozamiento
C     dUx(i,j,2) -> dUx en nodo i, componente j, corregido con rozamiento 

C     Other variables
      double precision g                    ! gravity acceleration
      double precision cfl                  ! Couriant-Friedichs-Levy number
      double precision k                    ! earth pressure coefficient
      double precision dt, t, t_end         ! time step, current time, simulation time
      double precision t_w, t_wn            ! 'writting' period, next 'writting' time

      double precision h_min
      double precision suma0,suma
      integer flag,l_print
      double precision fric1,fric2, fr_hmin, fr_qmin
      common /friction/ fric1, fric2,fr_hmin,fr_qmin
      common /genvars/ g,k,cfl,h_min,l_print
      character*40 achar
      integer c_sup,c_inf,c_izda,c_dcha,ncall
      double precision thrsld,mu,Ek,Fr,mumax
      integer fin
      integer ncp ! Number of cells Z > 0
      
      fin = 0
      h_min = 1.d-2
      nvar = 3
      k = 1.d0

      open(1,file='params.txt')
      read(1,*) fric1
      read(1,*) fric2
C      read(1,*) k
      close(1)

      fr_hmin = 1.d-2
      fr_qmin = 1.d-8

      flag = 0
      t = 0.d0
      t_end = 300.0
      t_w = 5.d0
      t_wn = t_w
      l_print = 50
      
      cfl = 0.9
      g = 9.8


C============================C
C                            C
C     Allocate arrays        C
C                            C
C============================C

      open(unit=1,file='malla0.txt')
      read(1,*) achar,nvert
      read(1,*) achar, ncx
      read(1,*) achar, ncy
      close(1)
      nnx = ncx + 1
      nny = ncy + 1

      nc = (nnx-1)*(nny -1)
      nn = nnx * nny
      nvw = nnx*(nny-1)
      nhw = (nnx-1)*nny
      nw = nvw + nhw

      allocate(U(nc,nvar),F1(nc,nvar),G1(nc,nvar))
      allocate(Z(nc),angle(nc,2))


      allocate(dUn(nw,nvar,2),dUp(nw,nvar,2))
      allocate(node(nn,2),celda(nc,2))
      allocate(hw(nhw,2),vw(nvw,2))
      allocate(hwn(nhw,2),vwn(nvw,2))
      allocate(cellw(nc,nvert),cellnod(nc,nvert))


C=============================================C
C                                             C
C     Read mesh and initial conditions        C
C                                             C
C=============================================C
 
      
      call read_mesh2D(node,celda,hw,vw,hwn,vwn,cellw,cellnod,Z,U,ncp)


CCCCCCCCCCCCCCCCCCCCCCCCC
C                       C
C     Compute angles    C
C                       C
CCCCCCCCCCCCCCCCCCCCCCCCC
    
      do i = 1, nc
         
         c_inf = hw(cellw(i,1),1)
         c_sup = hw(cellw(i,2),2)
         c_izda =vw(cellw(i,4),1)
         c_dcha =vw(cellw(i,3),2)

         if (c_izda.eq.0) then
            angle(i,1) = atan((Z(c_dcha)-Z(i))/
     +           (celda(c_dcha,1)-celda(i,1)))
         else if (c_dcha.eq.0) then
            angle(i,1) = atan((Z(i)-Z(c_izda))/
     +           (celda(i,1)-celda(c_izda,1)))
         else
            angle(i,1) = atan((Z(c_dcha)-Z(i))/
     +           (celda(c_dcha,1)-celda(i,1)))
            angle(i,1) = (angle(i,1) +  atan((Z(i)-Z(c_izda))/
     +           (celda(i,1)-celda(c_izda,1))))*0.5d0
         endif

         if (c_inf.eq.0) then
            angle(i,2) = atan((Z(c_sup)-Z(i))/
     +           (celda(c_sup,2)-celda(i,2)))
         else if (c_sup.eq.0) then
            angle(i,2) = atan((Z(i)-Z(c_inf))/
     +           (celda(i,2)-celda(c_inf,2)))
         else
            angle(i,2) = atan((Z(c_sup)-Z(i))/
     +           (celda(c_sup,2)-celda(i,2)))
            angle(i,2) = (angle(i,2) + atan((Z(i)-Z(c_inf))/
     +           (celda(i,2)-celda(c_inf,2))))*0.5d0
         endif

      enddo

C===========================================C
C                                           C
C     Go to conservative formulation        C
C                                           C
C===========================================C

      do j = 1, nc
         U(j,1) = U(j,1) * 1.15d0
      enddo

      suma0 = 0.d0
      do j = 1, nc
         do i = 2, nvar
            U(j,i) = U(j,i)*U(j,1)
         enddo
         do i = 2, nvar
            if (U(j,1).gt.h_min) then
               F1(j,i) = U(j,i)*U(j,2)/U(j,1)
            else

               F1(j,i) = 0.d0
            endif
            if (U(j,1).gt.h_min) then
               G1(j,i) = U(j,i)*U(j,3)/U(j,1)
            else

               G1(j,i) = 0.d0
            endif
            
         enddo
         suma0 = suma0 + U(j,1)
         F1(j,1) = U(j,2)
         G1(j,1) = U(j,3)
      enddo

C!C      open(1000,file='out')
C!C      open(33,file='Ekin')
C===========================================C
C                                           C
C                 Time Loop                 C
C                                           C
C===========================================C
      l = 0
      ncall = 0
!      call write_plt(U,Z,celda,ncall,fin,ncp)
      ncall = ncall + 1
C      stop

      do while (t < t_end)
         if (mod(l,l_print).eq.0) then
!            write(*,*) '===================================='
!            write(*,*) 'paso l=',l
         endif
         if (abs(t-t_wn).lt.1.d-5) then
!            call write_plt(U,Z,celda,ncall,fin,ncp)
            ncall = ncall + 1
            t_wn = t_wn + t_w
         endif
         

C         if (l.gt.36) stop

         l= l +1

         call calc_paredes2D(celda,node,hw,vw,hwn,vwn,cellw,U,F1
     +        ,G1,Z,angle,dUp,dUn,dt)
         
         if (t + dt.gt.t_wn) dt = t_wn - t
C         if (t + dt.gt.t_end) dt = t_end - dt

         call method2D(node,hw,vw,hwn,vwn,cellw,U
     +     ,dUp,dUn,dt,l)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Comprobar obstáculos    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                            C
C     Si la celda X comparte pared k         C
C     con una seca Y con z_Y > z_X + h_X     C
C     entonces el flujo en la direción       C
C     normal a k es cero                     C
C                                            C
C                              C
C     Paredes horizontales     C
C                              C
         do i = 1, nhw
            lc = hw(i,1)        ! celda izda / inferior
            rc = hw(i,2)        ! celda dcha / superior
            thrsld = 1.d-6
C     if lc seca y rc mojada
            if ((U(lc,1).lt.thrsld).and.(U(rc,1).gt.thrsld)) then
               if (Z(rc)+U(rc,1).lt.Z(lc)) U(rc,3) = 0.d0
            endif
C     if lc mojada y rc seca
            if ((U(lc,1).gt.thrsld).and.(U(rc,1).lt.thrsld)) then
               if (Z(lc) + U(lc,1).lt.Z(rc)) U(lc,3) = 0.d0
            endif
         enddo
C                            C
C     Paredes verticales     C
C                            C
         do i = 1, nvw
            lc = vw(i,1)        ! celda izda / inferior
            rc = vw(i,2)        ! celda dcha / superior
            thrsld = 1.d-6
C     if lc seca y rc mojada
            if ((U(lc,1).lt.thrsld).and.(U(rc,1).gt.thrsld)) then
               if (Z(rc)+U(rc,1).lt.Z(lc)) U(rc,2) = 0.d0
            endif
C     if lc mojada y rc seca
            if ((U(lc,1).gt.thrsld).and.(U(rc,1).lt.thrsld)) then
               if (Z(lc) + U(lc,1).lt.Z(rc)) U(lc,2) = 0.d0
            endif
         enddo

         t = t + dt

         suma = 0.d0
         Ek = 0.d0
         mumax = 0.d0
         flag = 1
         do j = 1, nc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SALIDA DEL ÁREA CENTRAL => FIN DE LA SIMULACION    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            c_inf = hw(cellw(j,1),1)
            c_sup = hw(cellw(j,2),2)
            c_izda =vw(cellw(j,4),1)
            c_dcha =vw(cellw(j,3),2)
            flag = c_inf*c_sup*c_izda*c_dcha
!!C            if ((U(j,1).gt.1.d-10).and.(flag.eq.0)) then
!!!C               write(*,*) 'Se derramó fuera'
!!!C               write(*,*) 'E = -100000'
!!C               open(111,file='Nash_Sutcliffe',POSITION='APPEND')
!!C               write(111,*) 1000000.0,fric1, fric2
!!C               close(111)
!!!C               write(*,*) 1000000.0
!!C               stop
!!C            endif

CcccccccccccccccccccccccccccccccccccccccccccccccccC
C     CALADOS NEGATIVOS : SALIDA POR PANTALLA     C
CcccccccccccccccccccccccccccccccccccccccccccccccccC

CC            if (U(j,1).lt.0.d0) then
!CC               write(*,*) 'Celda',j,'(',celda(j,1:2),'):',U(j,1)
CC            endif

CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMC
C     CORRECCIÓN: CALADOS MUY PEQUEÑOS => CAUDALES DETENIDOS    C
CWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWC
            if (U(j,1).lt.h_min) then
               U(j,2:3) = 0.d0
            else
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Energía cinética     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
               mu = (U(j,2)*U(j,2) + U(j,3)*U(j,3))/(U(j,1)*U(j,1))
               if (mu.gt.mumax) mumax = mu
               Ek = Ek + 0.5 * mu * U(j,1)
            endif
            do i = 2, nvar
               if (U(j,1).gt.h_min) then
                  F1(j,i) = U(j,i)*U(j,2)/U(j,1)
               else
                  F1(j,i) = 0.d0
               endif
               if (U(j,1).gt.h_min) then
                  G1(j,i) = U(j,i)*U(j,3)/U(j,1)
               else
                  G1(j,i) = 0.d0
               endif
              
            enddo
            if (U(j,1).gt.h_min) then
               mu = sqrt(mu)
               Fr = mu/sqrt(g*U(j,1))
            endif
            suma = suma + U(j,1)
            F1(j,1) = U(j,2)
            G1(j,1) = U(j,3)
         enddo

C!C         write(33,*) t,l,Ek,mumax
C!C         write(100,*) t,'Error de masa =',(suma-suma0)/suma0*100.d0,'%'
C!C         write(100,*) '          ',suma
         if (mod(l,l_print).eq.0) then
!            write(*,*) l,t,dt
!            write(*,*) '     ',Ek,(suma-suma0)/suma0*100.d0,'%'
         endif
         if ((t.gt.200).and.(Ek.lt.100.d0)) then
!            write(*,*) 'Stopped mass'
            fin = 1
!            call write_plt(U,Z,celda,ncall,fin,ncp)
!            call comparador(U,celda,suma)
            goto 219
         endif
      enddo
      fin = 1
!      call write_plt(U,Z,celda,ncall,fin,ncp)
!      call comparador(U,celda,suma)
      
 219  continue
C!C      close(33)
C!C      close(1000)

C!C      open(33,file='final')
C!C      do i = 1, nc
C!C         write(33,*) U(i,1)
C!C      enddo
C!C      close(33)
      deallocate(U)
      deallocate(F1,G1)
      deallocate(Z,angle)


      deallocate(dUp)
      deallocate(dUn)
      deallocate(node,celda)
      deallocate(hw,vw)
      deallocate(hwn,vwn)
      deallocate(cellw)
      
!!!      write(*,*) 'hecho'

      end
