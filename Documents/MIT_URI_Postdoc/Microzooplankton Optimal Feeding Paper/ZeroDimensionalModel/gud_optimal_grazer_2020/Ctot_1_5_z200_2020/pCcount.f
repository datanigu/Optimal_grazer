c ****** Calculating carbon content per cell and size class and number of cells per size class

         subroutine pCcount(pnum,pr,pv,pCcoeff,pCexp,
     &       pCsize,pCcell,pcount)

c *** Parameters for carbon concentration

        implicit none
        integer pnum
        real*8 pr(pnum) !cell radius
        real*8 pv(pnum)! cell volumne, m^3
        real*8 pCcoeff! coefficient for carbon per size class, g/m^3
        real*8 pCexp! dimensionless, exponent for carbon per size class
        real*8 pCsize(pnum)! carbon concentration in each size class, gC/m^3        
        real*8 pCcell(pnum)! carbon per cell, gC/cell
        real*8 pcount(pnum)! number of cells per vol per size class, cells/m^3
        integer i


c *** Calculating carbon per cell, carbon per size class, cells per size class
      do 50 i = 1,pnum
      pCcell(i) = (0.216*((pv(i)*10.0**18)**0.939))*10.0**(-12)
      pCsize(i) = pCcoeff*(pr(i)**pCexp)
      pcount(i) = pCsize(i)/pCcell(i)
50     continue

c       print*,'pCcell=',pCcell
c       print*,'pCsize=',pCsize
c       print*,'pcount=',pcount


       return
       end


