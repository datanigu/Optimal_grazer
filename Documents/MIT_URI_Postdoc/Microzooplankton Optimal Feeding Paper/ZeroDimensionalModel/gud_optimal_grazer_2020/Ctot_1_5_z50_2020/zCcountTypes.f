c **************** Calculating carbon and abundance related grazer variables *****
c *** Use this subroutine to calculate the amount of carbon per cell and per size class
c *** as well as the number of cells per size class


       subroutine zCcountTypes(fnum,zmax,zrmx,zvmx,zCcoeff,zCexp,
     &      functype,
     &     zCcell,zCsize,zcount)

c *** Declaring variables ***
      implicit none
      integer fnum ! number of functional groups
      integer zmax !absolut max number of size classes
      real*8 zvmx(zmax,fnum) !biovolume, m^3
      real*8 zrmx(zmax,fnum) !cell radius, m
      real*8 zCcoeff(fnum)! coefficient for carbon within each size class
      real*8 zCexp(fnum) ! exponent for carbon per size class
      integer functype(fnum)! type of functional group 
c               (1=dino, 2= ciliate, 3=generic microzoopl)
      real*8 zCsize(zmax,fnum)! carbon per size class, gC/m^3/size class
      real*8 zCcell(zmax,fnum) !carbon per cell, gC/cell
      real*8 zcount(zmax,fnum) !cells per size class, cells/zie class
      integer l, j1,j2 !indices for loops

c *** Calculating carbon-associated values
    
      do 350 l = 1,fnum
        if (functype(l) .eq. 1) then !dinoflagellate 
         do 375 j1 =1,zmax
      zCcell(j1,l) = (0.4436*((zvmx(j1,l)*10.0**18.0)**0.864))
     &              *(10.0**(-12.0))!(based on Menden-Deuer and Lessard 2000)
      zCsize(j1,l) = zCcoeff(l)*(zrmx(j1,l)**zCexp(l))! carbon concentration in each size class, gC/m^3
c         zcount(j1,l) = zCsize(j1,l)/zCcell(j1,l)
375        continue




       elseif (functype(l) .eq. 2 .or. functype(l) .eq. 3) then !ciliate or generic microzoopl. Yes, same equation
        do 450 j2 = 1,zmax
      zCcell(j2,l)=(0.216*((zvmx(j2,l)*10.0**18.0)**0.939))
     &          *(10.0**(-12.0))! carbon per cell, gC/cell
      zCsize(j2,l) = zCcoeff(l)*(zrmx(j2,l)**zCexp(l))! carbon concentration in each size class, gC/m^3
c        zcount(j2,l) = zCsize(j2,l)/zCcell(j2,l)
450     continue
        endif !for functype

350   continue

c      print*,'carbon per grazer cell=',zCcell
c      print*,'carbon conc. within each size class (zCsize)=',zCsize
ccc      print*,'cell count per size class=',zcount
c
c
      return
      end
 



