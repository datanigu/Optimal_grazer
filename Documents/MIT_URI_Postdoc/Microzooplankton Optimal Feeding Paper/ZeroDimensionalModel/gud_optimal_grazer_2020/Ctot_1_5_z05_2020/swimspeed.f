c        Swim Speed Changes with Grazing Intake

       subroutine swimspeed(fnum,zmax,maxind,minind,
     &     vmax,vmin, R0, Rv, g,v)

         implicit none
         integer fnum ! number of functional groups
         integer zmax ! max number of size classes
         real*8 maxind(fnum) !index for max size class
         real*8 minind(fnum) ! indiex for min size class
         real*8 vmax(zmax, fnum) !max allowed swimming speed, m/s
         real*8 vmin(zmax, fnum) ! min allowed swimming speed, m/s
         real*8 R0(zmax,fnum)  !basal respiration, 1/s
         real*8 Rv(zmax,fnum)  ! respiration from swimming, 1/s
         real*8 g(zmax,fnum)  ! grazing rate at a specific time point, 1/s
         real*8 v(zmax,fnum)  ! resulting swim speed, m/s
         integer i, f
         integer ind1, ind2
     

        do 150 f = 1,fnum
           ind1 = minind(f)
           ind2 = maxind(f)
         do 100 i = ind1, ind2
          if(g(i,f) .le. R0(i,f)) then
            v(i,f) = vmin(i,f)
         else if (g(i,f).gt.R0(i,f).and.
     &          g(i,f) .lt. R0(i,f)+Rv(i,f)) then
            v(i,f) = (vmax(i,f) - vmin(i,f))/Rv(i,f)*(g(i,f)-R0(i,f)) 
     &              + vmin(i,f);
         else if (g(i,f) .ge. (R0(i,f)+Rv(i,f))) then
            v(i,f) = vmax(i,f)
         endif
100      continue
150      continue

c       print*,'v IN subroutine',v

       end! subroutine

