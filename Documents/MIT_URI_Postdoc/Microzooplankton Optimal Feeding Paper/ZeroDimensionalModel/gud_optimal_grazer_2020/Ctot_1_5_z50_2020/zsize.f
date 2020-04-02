c ******** Grazer Sizes *************
c *** Use this subroutine to calculate the sizes (radius and volume) 
c *** of grazers for each functional group
c ** The output will be in matrix form, with each row being the max number 
c ** of size classes, and each column being a different functional group


      subroutine zsize(zmax,fnum,pi,sizemax,sizemin,zsizewidth,
     &   zrmx,znumf,zvmx,minind, maxind)

c *** Declaring variables of interest
      implicit none
      integer zmax !absolute max number of size classes to have (just a placeholder)
      integer fnum
      double precision pi
      real*8 sizemin(fnum)
      real*8 sizemax(fnum)
      real*8 zsizewidth ! log base size width, usually 2 or 10
      real*8 maxexp !exponent to find grazer size classes
      integer maxexpint!rounded integer form of exponent for do loop
      real*8 zrall(zmax)! all of the potential size classes
      real*8 zrmx(zmax,fnum)
      real*8 minind(fnum)!index of min size class out of zrall
      real*8 maxind(fnum)!index of max size class out of zrall
      integer znum ! max number of size classes out of all func gps
      integer znumf(fnum) !vector of number of size classes for each functional group
      real*8 zvmx(zmax,fnum)!grazer volume, m^3
      integer ind1, ind2 !for use in indexing WITHIN loop
      integer i, j, k !loop indices
 

c ************ Calculating all size classes across all functional groups
      maxexp = maxval(sizemax)/minval(sizemin)
      maxexp = log(maxexp)/log(zsizewidth)!find exponent for exopential size classes
      maxexpint = nint(maxexp)! round to nearest integer
     
      do 100 i = 1,maxexpint!+1
      zrall(i) = minval(sizemin)*(zsizewidth**(i-1))!all possible grazer radii
100   continue

c      print*, maxexpint, maxexp

c Finding size range for each functional group
       do 150 j = 1,fnum
       minind(j) = minloc(abs(zrall-sizemin(j)),1)
       maxind(j) = minloc(abs(zrall-sizemax(j)),1)
150    continue


c Finding max number of size classes out of all functional groups and putting
c those sizes in a matrix
       znumf = maxind - minind
       znumf = znumf +1
c       print*, 'znumf = ',znumf      
       zrmx = 0.d0
       do 300 k = 1,fnum
            ind1 = minind(k)
            ind2 = maxind(k)
            zrmx(ind1:ind2,k) = zrall(ind1:ind2)
            zvmx(:,k) =4.0/3.0*pi*(zrmx(:,k)**3)!volume, m^3

300    continue



c      print*, 'zrall=',zrall
c      print*, 'number size classes in each func group=',znumf
c      print*, 'minind=',minind
c      print*, 'maxind=',maxind
c      print*, 'zrmx col 1 = ', zrmx(:,1)
c      print*, 'zrmx col 2 = ', zrmx(:,2)
c       print*,'zrmx', zrmx
c      print*, 'zvmx col 1=',zvmx(:,1)
c      print*, 'zvmx col 2=',zvmx(:,2)
c
c
c

      return
      end
