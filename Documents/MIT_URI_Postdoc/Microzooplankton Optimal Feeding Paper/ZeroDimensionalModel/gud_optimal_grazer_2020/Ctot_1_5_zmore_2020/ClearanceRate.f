c **************** Calculating clearance rate *****************

       subroutine ClearanceRate(znum, fnum,pnum,minind,maxind,
     &   gmax, zCcell,beta, h,
     &  clearance)

          implicit none
          integer fnum! number of functional groups
          integer znum ! max number of grazer size classes
          integer pnum ! number of phytoplankton size classes
          real*8 minind(fnum) ! index for first size class for each func gp
          real*8 maxind(fnum) ! index for last size class for each func gp
          real*8 gmax(znum,fnum) ! max grazing rate
          real*8 zCcell(znum,fnum) ! carbon in each grazer size class
          real*8 beta(znum,pnum,fnum) ! encounter kernel for each pred-prey combo
          real*8 h(znum,pnum,fnum) ! handling time for each pred-prey combo
          real*8 clearance(znum,pnum,fnum) ! clearance rate
          integer i ! indices for loops
          integer j
          integer k
          integer i1, i2
             
          do 100 i = 1,fnum
             i1 = minind(i)
             i2 = maxind(i)
        
c             print*,'minind(i)',minind(i)
c             print*,'maxind(i)',maxind(i)
c             print*,'i1',i1
c             print*,'i2',i2
c             print*,'gmax(i1:i2,i)',gmax(i1:i2,i)
c             print*,'zCcell(i1:i2,i)',zCcell(i1:i2,i)
c             print*,'beta(i1:i2,j,i)',beta(i1:i2,j,i)
c             print*,'h(i1:i2,j,i)',h(i1:i2,j,i)
          do 300 k = i1,i2
            do 200 j = 1,pnum
          clearance(k,j,i) = 
     &         gmax(k,i)
     &         /(zCcell(k,i)/beta(k,j,i)/h(k,j,i))
          
200         continue
300        continue
100       continue

       end !subroutine
