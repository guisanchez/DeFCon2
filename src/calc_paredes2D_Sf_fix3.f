
      subroutine calc_paredes2D(celda,node,hw,vw,hwn,vwn,cellw,U,F1
     +     ,G1,Z,angle,dUp,dUn,dt)

      implicit none

      integer nc,nn,nnx,nny,nvert,nvar
      integer nvw,nhw,nw,l_print
      double precision g,k,cfl,h_min
      common /n_cell/ nc,nn,nnx,nny,nvert,nvar
      common /n_wall/ nvw,nhw
      common /genvars/ g,k,cfl,h_min,l_print

      
      integer i,j,l
      double precision node(nn,2),celda(nc,2)
      integer cellw(nc,nvert)
      double precision U(nc,nvar),F1(nc,nvar),G1(nc,nvar)

      double precision Z(nc),dz
      double precision dUp(nhw+nvw,nvar,2),dUn(nhw+nvw,nvar,2)
      double precision angle(nc,2),theta(2),c2t
      double precision l_1m,l_2m,l_3m

      double precision Pi(3,3),P_invi(3,3),Opi(3,3),Omi(3,3),MP(3,3)
      double precision MPn(nhw+nvw,3,3), MPp(nhw+nvw,3,3)
      integer hw(nhw,2),vw(nvw,2),hwn(nhw,2),vwn(nvw,2)
      integer pared(nhw+nvw,2)
      double precision pml(nhw+nvw)  ! pared:máximo_lambda
      double precision ubar,vbar,cbar,hbar
      double precision lambdamax,dt,dtmin
      double precision dh,dF(nvar),dG(nvar),S(nvar),u1,u2,u3
      integer rc,lc

      double precision qps     ! q-predicted-sign (to evaluate friction sign)
     
      double precision Ai,dx,dy,dL
      integer pin,ps,pli,pld,sc
      
      double precision Sfric(nhw+nvw),fac(nhw+nvw),fac0
      double precision Sfmax(nhw+nvw),ex_n,ex_p,suma
      double precision fwall(nhw+nvw)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Corrección "entrópica"     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision visc, l1_l,l1_r,l2_r,l2_l,l3_l,l3_r,vaux
      integer nup

      double precision hup


      sc = 0

      nw = nhw + nvw ! Num. walls = Num. horizontal walls + Num. vertical walls
      pared(1:nhw,1:2) = hw(1:nhw,1:2)
      pared(1+nhw:nw,1:2) = vw(1:nvw,1:2)

      
      lambdamax = 0.0
      dtmin = 1000.0


      nup = 0
      hup = 1.d-6

C=============================C
C                             C
C     CALCULO POR PAREDES     C
C                             C
C=============================C

      rc = 0
      lc = 0

      do i = 1, nw

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        Si no estamos en el borde         C  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         l1_l = 0.d0
         l2_l = 0.d0
         l3_l = 0.d0
         l1_r = 0.d0
         l2_r = 0.d0
         l3_r = 0.d0
         visc = 0.d0

         if (pared(i,1)*pared(i,2).ne.0) then !En el borde pared(i,1 ó 2)=0 
            
            lc = pared(i,1)
            rc = pared(i,2)

            ubar = 0.d0
            vbar = 0.d0
            theta(1) = (angle(rc,1) + angle(lc,1)) / 2.0
            theta(2) = (angle(rc,2) + angle(lc,2)) / 2.0
            c2t = 1.d0/(1.d0 + tan(theta(1))*tan(theta(1)) +
     +           tan(theta(2))*tan(theta(2)))
C            c2t = cos(theta(1))*cos(theta(1)) +
C     +           cos(theta(2))*cos(theta(2))
            hbar = (U(lc,1) + U(rc,1)) / 2.0
            cbar = sqrt(hbar*g*k*c2t)
            if ((U(lc,1).ge.1.d-6).and.(U(rc,1).ge.1.d-6)) then
               ubar = U(lc,2)/sqrt(U(lc,1)) + U(rc,2)/sqrt(U(rc,1))
               ubar = ubar /(sqrt(U(lc,1)) + sqrt(U(rc,1)))
               vbar = U(lc,3)/sqrt(U(lc,1)) + U(rc,3)/sqrt(U(rc,1))
               vbar = vbar /(sqrt(U(lc,1)) + sqrt(U(rc,1)))
            else if ((U(lc,1).lt.1.d-6).and.(U(rc,1).ge.1.d-6)) then
               ubar = U(rc,2)/sqrt(U(rc,1))
               ubar = ubar /(sqrt(U(rc,1)))
               vbar = U(rc,3)/sqrt(U(rc,1))
               vbar = vbar /(sqrt(U(rc,1)))
            else if ((U(lc,1).ge.1.d-6).and.(U(rc,1).lt.1.d-6)) then
               ubar = U(lc,2)/sqrt(U(lc,1)) 
               ubar = ubar /(sqrt(U(lc,1))) 
               vbar = U(lc,3)/sqrt(U(lc,1)) 
               vbar = vbar /(sqrt(U(lc,1))) 
            else                
               ubar = 0.d0
               vbar = 0.d0
            endif
         else ! i.e., si la pared i está en el borde (la celda 0 es vecina)
            ubar = 0.d0
            vbar = 0.d0
            cbar = 0.d0
            
         endif


         if (i.le.nhw) then     ! Si la pared es horizontal
            l_1m = vbar - cbar 
            l_2m = vbar + cbar 
            l_3m = vbar             
         else                   ! Si la pared es vertical
            l_1m = ubar - cbar 
            l_2m = ubar + cbar 
            l_3m = ubar 
         endif
         
         fwall(i) = max(abs(l_1m),abs(l_2m),abs(l_3m))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     If fwall(i) < fwall_min -> no flux across wall 'i'    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
         
         vaux = 0.d0
         
         if (pared(i,1)*pared(i,2).ne.0) then

            if (i.le.nhw) then
               vaux = U(lc,3)
            else
               vaux = U(lc,2)
            endif
            
            if (U(lc,1).gt.1.d-6) then
               l1_l = vaux/U(lc,1) -
     +              sqrt(((cos(angle(lc,1))*cos(angle(lc,1)) +
     +              cos(angle(lc,2)))*cos(angle(lc,2)))*g*k*U(lc,1))
               l2_l = vaux/U(lc,1) +
     +              sqrt(((cos(angle(lc,1))*cos(angle(lc,1)) +
     +              cos(angle(lc,2)))*cos(angle(lc,2)))*g*k*U(lc,1))
               l3_l = vaux/U(lc,1) 
            endif

            if (i.le.nhw) then
               vaux = U(rc,3)
            else
               vaux = U(rc,2)
            endif
            if (abs(U(rc,1)).gt.1.d-6) then
               l1_r = vaux/U(rc,1) -
     +              sqrt(((cos(angle(rc,1))*cos(angle(rc,1)) +
     +              cos(angle(rc,2)))*cos(angle(rc,2)))*g*k*U(rc,1))
               l2_r = vaux/U(rc,1) +
     +              sqrt(((cos(angle(rc,1))*cos(angle(rc,1)) +
     +              cos(angle(rc,2)))*cos(angle(rc,2)))*g*k*U(rc,1))
               l3_r = vaux/U(rc,1)
            endif
         
         endif
         
         
         if ((l1_l.lt.0.d0).and.(l1_r.gt.0.d0)) 
     +        visc = l1_r - l1_l
         if ((l2_l.lt.0.d0).and.(l2_r.gt.0.d0)) 
     +        visc = max(visc, l2_r - l2_l)
         if ((l3_l.lt.0.d0).and.(l3_r.gt.0.d0))
     +        visc = max(visc, l3_r - l3_l)

         pml(i) = max(abs(l_1m),abs(l_2m),abs(l_3m))
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Definiendo matrices P y P^{-1}    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         if (i.le.nhw) then
CCCCCCCCCCCCCCCCCCCCCCCCCC
C  Paredes horizontales  C     
CCCCCCCCCCCCCCCCCCCCCCCCCC     
            
            if (abs(cbar).gt.1.d-11) then
               Pi(1,1:2) = 1.d0
               Pi(1,3) = 0.d0
               Pi(2,1:2) = ubar
               Pi(2,3) = 1.d0
               Pi(3,1) = l_1m
               Pi(3,2) = l_2m
               Pi(3,3) = 0.d0

               P_invi(1,1) = l_2m/(2.d0*cbar)
               P_invi(1,2) = 0.d0
               P_invi(1,3) = -1.d0/(2.0*cbar)
               P_invi(2,1) = -l_1m/(2.d0*cbar)
               P_invi(2,2) = 0.d0
               P_invi(2,3) = 1.d0/(2.d0*cbar)
               P_invi(3,1) = -ubar
               P_invi(3,2) = 1.d0
               P_invi(3,3) = 0.d0
            else
               Pi(1,1) = 1.d0
               Pi(1,2:3) = 0.d0
               Pi(2,1) = 0.d0
               Pi(2,2) = 1.d0
               Pi(2,3) = 0.d0
               Pi(3,1:2) = 0.d0
               Pi(3,3) = 1.d0

               P_invi(1,1) = 1.d0
               P_invi(1,2:3) = 0.d0
               P_invi(2,1) = 0.d0
               P_invi(2,2) = 1.d0
               P_invi(2,3) = 0.d0
               P_invi(3,1:2) = 0.d0
               P_invi(3,3) = 1.d0
            endif
            
         else
CCCCCCCCCCCCCCCCCCCCCCCC
C  Paredes verticales  C
CCCCCCCCCCCCCCCCCCCCCCCC
            
            if (abs(cbar).gt.1.d-11) then
               Pi(1,1:2) = 1.d0
               Pi(1,3) = 0.d0
               Pi(2,1) = l_1m
               Pi(2,2) = l_2m
               Pi(2,3) = 0.d0
               Pi(3,1:2) = vbar
               Pi(3,3) = 1.d0  

               P_invi(1,1) = l_2m/(2.d0*cbar)
               P_invi(1,2) = -1.d0/(2.0*cbar)
               P_invi(1,3) = 0.d0
               P_invi(2,1) = -l_1m/(2.d0*cbar)
               P_invi(2,2) = 1.d0/(2.d0*cbar)
               P_invi(2,3) = 0.d0
               P_invi(3,1) = -vbar
               P_invi(3,2) = 0.d0
               P_invi(3,3) = 1.d0
            else
               Pi(1,1) = 1.d0
               Pi(1,2:3) = 0.d0
               Pi(2,1) = 0.d0
               Pi(2,2) = 1.d0
               Pi(2,3) = 0.d0
               Pi(3,1:2) = 0.d0
               Pi(3,3) = 1.d0

               P_invi(1,1) = 1.d0
               P_invi(1,2:3) = 0.d0
               P_invi(2,1) = 0.d0
               P_invi(2,2) = 1.d0
               P_invi(2,3) = 0.d0
               P_invi(3,1:2) = 0.d0
               P_invi(3,3) = 1.d0
            endif            
         endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Definiendo matrices O^{+} y O^{-}     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
           
         if (l_1m.ne.0.) then
            Opi(1,1) = 0.5d0 * (1.d0 + sign(1.d+0,l_1m))
            Omi(1,1) = 0.5d0 * (1.d0 - sign(1.d+0,l_1m))
         else
            Opi(1,1) = 0.5d0
            Omi(1,1) = 0.5d0
         endif
         Opi(1,2:3) = 0.d0
         Opi(2,1) = 0.d0
         Omi(1,2:3) = 0.d0
         Omi(2,1) = 0.d0
         if (l_2m.ne.0.) then
            Opi(2,2) = 0.5d0 * (1.d0 + sign(1.d+0,l_2m))
            Omi(2,2) = 0.5d0 * (1.d0 - sign(1.d+0,l_2m))
         else
            Opi(2,2) = 0.5d0
            Omi(2,2) = 0.5d0
         endif
         Opi(2,3) = 0.d0
         Omi(2,3) = 0.d0
         Opi(3,1:2) = 0.d0
         Omi(3,1:2) = 0.d0
         if (l_3m.ne.0.d0) then
            Opi(3,3) = 0.5d0 * (1.d0 + sign(1.d+0,l_3m))
            Omi(3,3) = 0.5d0 * (1.d0 - sign(1.d+0,l_3m))
         else 
            Opi(3,3) = 0.5d0
            Omi(3,3) = 0.5d0
         endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Cálculo por paredes      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         if ((pared(i,1)*pared(i,2).ne.0).and.(fwall(i).gt.1.d-10)) then ! Si los dos vecinos existen (y hay sustancia)
            lc = pared(i,1)     ! Vecino 1
            rc = pared(i,2)     ! Vecino 2
            
            dh = U(lc,1) - U(rc,1) 
            if ((U(lc,1).lt.1.d-8).and.
     +           ((U(rc,1).lt.1.d-8))) then
               dz = 0.d0
            else
               dz = Z(lc) - Z(rc)
            endif
            
            u1 = 0.5d0*(U(lc,1) + U(rc,1)) 

CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXC
C           Calculos preliminares: Términos dU+           C
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Si la pared es horizontal:    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if (i.le.nhw) then

C     Longitud de la pared (que es horizontal)
               dL = abs(node(hwn(i,1),1)-node(hwn(i,2),1))
               
CCCCCCCCCCCCCCCCCCCCCCCCC
C     Término fuente    C
CCCCCCCCCCCCCCCCCCCCCCCCC               
               S(1) = 0.d0
               S(2) = 0.d0
               S(3) = g*c2t*u1*dL*(dh+dz)
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Diferencia flujos: (paredes horizontales)      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               
               dF(1:3) = 0.d0
               
               dG(1) = (G1(rc,1) - G1(lc,1))*dL
               dG(2) = (G1(rc,2) - G1(lc,2))*dL
               dG(3) = (G1(rc,3) - G1(lc,3))*dL
               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Si la pared es vertical:    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            else
               j = i - nhw
C     Longitud de la pared (que es vertical)
               dL = abs(node(vwn(j,1),2)-node(vwn(j,2),2))

CCCCCCCCCCCCCCCCCCCCCCCCC
C     Término fuente    C
CCCCCCCCCCCCCCCCCCCCCCCCC      

               S(1) = 0.d0
               S(2) = g*c2t*u1*dL*(dh+dz) 
               S(3) = 0.d0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Diferencia flujos: (paredes verticales)      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

               dF(1) = (F1(rc,1) - F1(lc,1))*dL
               dF(2) = (F1(rc,2) - F1(lc,2))*dL
               dF(3) = (F1(rc,3) - F1(lc,3))*dL

               dG(1:3) = 0.d0                  
            endif
!Dato para el cálculo de dt. dt = cfl*Ai/(max(l_m*dL)_horiz + max(l_m*dL)_vert)
            pml(i) = pml(i)*dL

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Independientemente de si la pared es | o -    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            call matprod(Pi,Opi,P_invi,MP,nvar)
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Predictor + a través de paredes     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            
            dUp(i,1:3,1:2) = 0.d0
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ajuste: sin fuente de presión desde celdas secas altas a     C
C                      celdas mojadas bajas                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

            if ((U(lc,1).lt.h_min).and.(Z(lc).gt.
     +           (Z(rc) + U(rc,1)))) then
               S(2:3) = 0.d0
               dF(2:3) = 0.d0
               dG(2:3) = 0.d0
            endif
            if ((U(rc,1).lt.h_min).and.(Z(rc).gt.
     +           (Z(lc) + U(lc,1)))) then
               S(2:3) = 0.d0  
               dF(2:3) = 0.d0
               dG(2:3) = 0.d0
            endif
            visc = visc*0.25d0*dL

            do j = 1, 3
               dUp(i,1,1) = dUp(i,1,1) + MP(1,j)*(S(j) - dF(j) - dG(j))
               dUp(i,2,1) = dUp(i,2,1) + MP(2,j)*(S(j) - dF(j) - dG(j))
               dUp(i,3,1) = dUp(i,3,1) + MP(3,j)*(S(j) - dF(j) - dG(j))
            enddo
            dUp(i,1:3,1) =dUp(i,1:3,1) - visc*(U(rc,1:3) - U(lc,1:3))
    

            if (i.le.nhw) then
               Sfmax(i) = S(3) - dF(3) - dG(3) 
            else
               Sfmax(i) = S(2) - dF(2) - dG(2) 
            endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Se almacena la matriz para su uso posterior    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            
            MPp(i,1:3,1:3) = MP(1:3,1:3)
            
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXC
C           Calculos preliminares: Términos dU-           C
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXC
C
C     Ya conozco (calculados previamente) los dF, dG y S, tanto 
C     si la pared es vertical como si es horizontal
C
            call matprod(Pi,Omi,P_invi,MP,nvar)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Predictor - a través de paredes     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC               

            dUn(i,1:3,1:2) = 0.d0

            do j = 1, 3
               dUn(i,1,1) = dUn(i,1,1) + MP(1,j)*(S(j) - dF(j) - dG(j))
               dUn(i,2,1) = dUn(i,2,1) + MP(2,j)*(S(j) - dF(j) - dG(j))
               dUn(i,3,1) = dUn(i,3,1) + MP(3,j)*(S(j) - dF(j) - dG(j))
            enddo
            dUn(i,1:3,1) = dUn(i,1:3,1) + visc*(U(rc,1:3) - U(lc,1:3))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Se almacena la matriz para su uso posterior    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            MPn(i,1:3,1:3) = MP(1:3,1:3)
         
C                                                             C
C     Corrección final: De donde no hay no se puede sacar     C
C     Si U(lc,1)<0.01 && dUp(i,1)<0 -> estamos sacando        C
C     sustancia de donde ya no queda: Paralizar los flujos    C
C                                                             C
C            if ((U(rc,1).lt.1.d-2).and.(dUp(i,1,1).lt.0.d0)) then
C               dUp(i,1:3,1) = 0.d0
C               dUn(i,1:3,1) = 0.d0
C               Sfmax(i) = 0.d0
CC               write(*,*) 'Sfmax modificado en pared',i
C            else if ((U(lc,1).lt.1.d-2).and.(dUn(i,1,1).lt.0.d0)) then
C               dUp(i,1:3,1) = 0.d0
C               dUn(i,1:3,1) = 0.d0
C               Sfmax(i) = 0.d0
C               write(*,*) 'Sfmax modificado en pared',i
C            endif

         else                   ! i.e., si no existen ambos vecinos
            dUn(i,1:nvar,1) = 0.d0 ! los flujos son inexistentes
            dUp(i,1:nvar,1) = 0.d0 ! (porque es pared frontera del sistema)
         endif
      enddo


C===================================C
C     Cálculo del paso temporal     C
C===================================C

      do i = 1, nc   ! Bucle en celdas

C     Pared inferior (horizontal)
         pin = cellw(i,1)
C     Nodos que delimitan pared pi: node(hwn(pi,1),1:2) <->node(hwn(pi,2),1:2)
         dx = abs(node(hwn(pin,1),1) - node(hwn(pin,2),1))
C     Pared lateral izquierda (vertical)
         pli = cellw(i,4)
C     Nodos que delimitan pared pl: node(vwn(pli,1),1:2) <->node(vwn(pli,2),1:2)
         dy = abs(node(vwn(pli,1),2) - node(vwn(pli,2),2))
C     Pared lateral derecha (vertical)
         pld = cellw(i,3)
C     Pared superior (horizontal) 
         ps = cellw(i,2)
     
         pli = pli + nhw
         pld = pld + nhw
C>>>>>>>>>>>>>><<<<<<<<<<<<<<C
C     Área de la celda i     C
C>>>>>>>>>>>>>><<<<<<<<<<<<<<C
         Ai = dx*dy
         dt = cfl*Ai/(max(pml(pli),pml(pld))+max(pml(pin),pml(ps)))
         if (dt.lt.dtmin) dtmin = dt
      enddo

      dt = dtmin

C/////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\C
C     Nuevo bucle: Cálculo de términos correctores (rozamiento)    C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////C


      do i = 1, nw
         if (pared(i,1)*pared(i,2).ne.0) then
            lc = pared(i,1)
            rc = pared(i,2)
            dh = U(lc,1) - U(rc,1) 
            if ((U(lc,1).lt.1.d-8).and.
     +           ((U(rc,1).lt.1.d-8))) then
               dz = 0.d0
            else
               dz = Z(lc) - Z(rc)
            endif
            u1 = 0.5d0*(U(lc,1) + U(rc,1)) 
            u2 = 0.5d0*(U(lc,2) + U(rc,2)) 
            u3 = 0.5d0*(U(lc,3) + U(rc,3)) 

C                               C
C     Cálculo del rozamiento    C
C                               C

            dx = (celda(lc,1)-celda(rc,1))*(celda(lc,1)-celda(rc,1))
            dx = dx +(celda(lc,2)-celda(rc,2))*(celda(lc,2)-celda(rc,2))
            dx = sqrt(dx)       ! Distancia entre centros de celdas
C     Longitud de la pared
C                          horizontal
            if (i.le.nhw) then
               dL = abs(node(hwn(i,1),1)-node(hwn(i,2),1))
            else
C                          vertical
               j = i- nhw
               dL = abs(node(vwn(j,1),2)-node(vwn(j,2),2))
            endif
            if ((U(lc,1).ge.1.d-6).or.(U(rc,1).ge.1.d-6)) then  
               call compute_Sf2D(u1,u2,u3,nvar,S,dx,dL,c2t,i,nhw)

C               if (((lc.eq.477).or.(rc.eq.477)).and.(i.gt.nhw)) then
                  
C                  write(*,*) 'call compute sf',S(1:3)
C               endif
               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CALCULANDO CORRECTOR   (+)  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@C
C     Dirección del flujo    C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@C 
CCCCCCCCCCCCCCCCCCCCCCCCCC
C  Paredes horizontales  C     
CCCCCCCCCCCCCCCCCCCCCCCCCC
               u2 = 0.5d0*(U(lc,2) + U(rc,2) + dt/(dL*dx)*(
     +                 dUp(cellw(lc,4)+nhw,2,1) + dUn(i,2,1) + 
     +                 dUp(i,2,1) + dUn(cellw(rc,3)+nhw,2,1)))
               u3 = 0.5d0*(U(lc,3) + U(rc,3) + dt/(dL*dx)*(
     +                 dUp(cellw(lc,1),3,1) + dUn(i,3,1) + 
     +                 dUp(i,3,1) + dUn(cellw(rc,2),3,1)))
               if (i.le.nhw) then
                  qps = u3
                  if (u2*u2+u3*u3.gt.1.d-20) then
                     S = -S*u3/sqrt(u2*u2+u3*u3)
                  else
                     S = 0.d0
                  endif
                  Sfmax(i) = Sfmax(i) + 0.5*dx*dL*(U(lc,3) + 
     +                 U(rc,3))/dt
               else
C                  j = i - nhw
                  qps = u2
                  if (u2*u2+u3*u3.gt.1.d-20) then
                  S = -S*u2/sqrt(u2*u2+u3*u3)
                  else
                     S = 0.d0
                  endif

                  Sfmax(i) = Sfmax(i) + 0.5*dx*dL*(U(lc,2) + 
     +                 U(rc,2))/dt
               endif

               if (abs(qps).gt.1.d-11) then
                  Sfmax(i) = sign(Sfmax(i),-qps)
                  do l = 1, nvar
                     if (S(l)*Sfmax(i).le.0.d0) then
                        S(l) = 0.d0
                     else if (S(l).gt.0.d0) then
                        S(l) = min(S(l),Sfmax(i))
                     else
                        S(l) = max(S(l),Sfmax(i))
                     endif
                  enddo
               else
                  S(1:nvar) = 0.d0
               endif
            else
               S(1:nvar) = 0.d0
            endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ajuste: sin fuente de presión desde celdas secas altas a     C
C                      celdas mojadas bajas                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
            if ((U(lc,1).lt.h_min).and.(Z(lc).gt.
     +           (Z(rc) + U(rc,1)))) S(1:3) = 0.d0   
            if ((U(rc,1).lt.h_min).and.(Z(rc).gt.
     +           (Z(lc) + U(lc,1)))) S(1:3) = 0.d0  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ajuste: Si no hay flujo a través de dicha pared,     C
C                no hay tampoco rozamiento                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C            if ((dUp(i,1,1).eq.0.d0).and.(dUp(i,2,1).eq.0.d0).and.
C     +           (dUp(i,3,1).eq.0.d0)) S(1:3) = 0.d0


            do j = 1, 3
               dUp(i,1,2) = dUp(i,1,2) + MPp(i,1,j)*S(j)
               dUp(i,2,2) = dUp(i,2,2) + MPp(i,2,j)*S(j)
               dUp(i,3,2) = dUp(i,3,2) + MPp(i,3,j)*S(j)
            enddo

C            if (((lc.eq.477).or.(rc.eq.477)).and.(i.gt.nhw)) then
C               write(*,*) 'wall',i,'cells',lc,rc, 'dUp',dUp(i,2,2)
C               write(*,*) 'Sf=',S(1:3)
C            endif

            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CALCULANDO CORRECTOR   (-)  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ajuste: sin fuente de presión desde celdas secas altas a     C
C                      celdas mojadas bajas                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
            if ((U(lc,1).lt.1.d-2).and.(Z(lc).gt.
     +           (Z(rc) + U(rc,1)))) S(1:3) = 0.d0   
            if ((U(rc,1).lt.1.d-2).and.(Z(rc).gt.
     +           (Z(lc) + U(lc,1)))) S(1:3) = 0.d0                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ajuste: Si no hay flujo a través de dicha pared,     C
C                no hay tampoco rozamiento                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C            if ((dUn(i,1,1).eq.0.d0).and.(dUn(i,2,1).eq.0.d0).and.
C     +           (dUn(i,3,1).eq.0.d0)) S(1:3) = 0.d0            

            do j = 1, 3
               dUn(i,1,2) = dUn(i,1,2) + MPn(i,1,j)*S(j)
               dUn(i,2,2) = dUn(i,2,2) + MPn(i,2,j)*S(j)
               dUn(i,3,2) = dUn(i,3,2) + MPn(i,3,j)*S(j)
            enddo

            ex_n = abs(dUn(i,1,2)) - abs(dUn(i,1,1))
            ex_p = abs(dUp(i,1,2)) - abs(dUp(i,1,1))
            
            if (max(ex_n,ex_p).gt.1.d-8) then

               ex_n = max(ex_n,ex_p)
               if (dUn(i,1,2).lt.0.d0) then
                  dUn(i,1,2) = dUn(i,1,2) + ex_n
                  dUp(i,1,2) = dUp(i,1,2) - ex_n
               else if (dUn(i,1,2).gt.0.d0) then
                  dUn(i,1,2) = dUn(i,1,2) - ex_n
                  dUp(i,1,2) = dUp(i,1,2) + ex_n
               endif
            endif

            if (i.le.nhw) then
               Sfric(i) = S(3)
            else
               Sfric(i) = S(2)
            endif
         else
            dUp(i,1:3,2) = 0.d0
            dUn(i,1:3,2) = 0.d0
         endif
      enddo

      goto 300

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Bucle en Paredes: Comprobar si corrección no supera a la predicción     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      fac(1:nw) = -1.d0

      do i = 1, nw
         if (pared(i,1)*pared(i,2).ne.0) then
C     En una pared tal que el flujo, calculado usando solo
C     la advección, tenga una dirección determinada, no puede
C     haber flujo en la dirección opuesta cuando se incluye
C     el rozamiento. Si |dUpred| < |dUcorr| -> recorregir dUcorr
C     fac = |dUpred|/|dUcorr|, 0 < fac < 1
C     dUcorr = MP * S: S -> S*fac
C     => dUcorr = MP * (S * fac)
            if (abs(dUp(i,1,2)).gt.abs(dUp(i,1,1))) then
               fac0 = abs(dUp(i,1,1))/abs(dUp(i,1,2))
            endif
            if (abs(dUn(i,1,2)).gt.abs(dUn(i,1,1))) then
               fac0 = min(abs(dUn(i,1,1))/abs(dUn(i,1,2)),fac0)
            endif

C     Se calcula el nuevo valor del rozamiento para esa pared

            if (i.le.nhw) then
               S(1) = 0.d0
               S(2) = 0.d0
               S(3) = Sfric(i)*fac0
            else
               S(1) = 0.d0
               S(2) = Sfric(i)*fac0
               S(3) = 0.d0
            endif
            dUp(i,1:3,2) = 0.d0
            do j = 1, 3
               dUp(i,1,2) = dUp(i,1,2) + MPp(i,1,j)*S(j)
               dUp(i,2,2) = dUp(i,2,2) + MPp(i,2,j)*S(j)
               dUp(i,3,2) = dUp(i,3,2) + MPp(i,3,j)*S(j)
            enddo
            dUn(i,1:3,2) = 0.d0
            do j = 1, 3
               dUn(i,1,2) = dUn(i,1,2) + MPn(i,1,j)*S(j)
               dUn(i,2,2) = dUn(i,2,2) + MPn(i,2,j)*S(j)
               dUn(i,3,2) = dUn(i,3,2) + MPn(i,3,j)*S(j)
            enddo
         else
            dUp(i,1:3,2) = 0.d0
            dUn(i,1:3,2) = 0.d0
         endif
      enddo

 300  continue
      goto 400
      suma = 0.d0
      do  i = 1, nw
         j = 1
         suma = suma + dUp(i,j,1) + dUn(i,j,1)
      enddo
      
      write(*,*) 'Sum. Predictors =',suma
      if (abs(suma).gt.1.d0) stop
      suma = 0.d0
      do  i = 1, nw
         
         j = 1
         suma = suma + dUp(i,j,2) + dUn(i,j,2)
         
      enddo
      write(*,*) 'Sum. Correctors =',suma
      
      
      if (sc.gt.0) then
         write(*,*) sc ,'paredes con flujo supercrítico'
      endif
 400  continue
      return 

      end
