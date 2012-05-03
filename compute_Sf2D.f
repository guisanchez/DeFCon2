      subroutine compute_Sf2D(h,qx,qy,nvar,Sf,dx,dL,c2t,w,nhw)

      implicit none

      integer w,nhw,nvar
      double precision Sff,Sf(nvar),c2t,dx,dL

      double precision c1,c2,u2,h,q,qx,qy,g,k,cfl
      double precision fric1,fric2,fr_hmin,fr_qmin
      common /friction/ fric1, fric2,fr_hmin,fr_qmin
      common /genvars/ g,k,cfl



      c1 = fric1
      c2 = fric2
      
      if (w.le.nhw) then 
         q = qy
      else
         q = qx
       endif

      if (h.gt.fr_hmin) then
         u2 = q*q/(h*h)
         Sff = c1 + c2*u2/h
      else
         u2 = 0.d0
         Sff = c1 !0.d0
      endif
      
      Sff = Sff*g*c2t*h*dx*dL*0.5d0

      Sf(1) = 0.d0

      if (w.le.nhw) then 
         Sf(2) = 0.d0
         Sf(3) = Sff
      else
         Sf(2) = Sff
         Sf(3) = 0.d0
      endif

C      if (abs(q).le.fr_qmin) then
CC      if (q.eq.0.d0) then
C         Sf = 0.d0
CC      else
CC         Sf = sign(Sf,q) !=Sf*q/abs(q)
C      endif

C      Sf(1:3) = 0.d0
      return 
      end
      
