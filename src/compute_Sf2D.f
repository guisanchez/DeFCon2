      subroutine compute_Sf2D(h,qx,qy,nvar,Sf,dx,dL,c2t,w,nhw)

      implicit none

      integer w,nhw,nvar
      double precision Sff,Sf(nvar),c2t,dx,dL
      integer frstyle
      double precision c1,c2,c3,u2,h,q,qx,qy,qm,g,k,cfl
      double precision fric1,fric2, fric3,fr_hmin,fr_qmin,rho
      common /friction/ fric1, fric2,fric3,frstyle,fr_hmin,fr_qmin,rho
      common /genvars/ g,k,cfl

      qm = sqrt(qx*qx+qy*qy)

      if (frstyle.le.2) then
         c1 = fric1
         c2 = fric2
      else 
         c1 = fric1
         c2 = fric2
         c3 = fric3
      endif

      if (w.le.nhw) then 
         q = qy
      else
         q = qx
      endif

C      if (qm.gt.1.d-6) then
C         c1 = c1 * q / qm
C      else
C         c1 = 0.d0
C      endif
      if (frstyle.eq.1) then
      
         if (h.gt.fr_hmin) then
            u2 = qm*qm/(h*h*h)
            Sff = c1 + c2*u2
         else
            u2 = 0.d0
            Sff = c1            !0.d0
         endif

      else if (frstyle.eq.2) then
         if (h.gt.fr_hmin) then
            u2 = 3.d0*c2*qm/(h*h)
            Sff = 1.d0/(h*rho)*(3.d0/2.d0*c1 + u2)
         else
            Sff = 0.d0
         endif
      else if (frstyle.eq.3) then
         if (h.gt.fr_hmin) then
            u2 = 3.d0*c3*qm/(h*h)
            Sff = c1 + 1.d0/(h*rho)*(3.d0/2.d0*c1 + u2)
         else
            Sff = c1
         endif
      endif

      Sff = Sff*g*c2t*h*dx*dL!*0.5d0

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
      
