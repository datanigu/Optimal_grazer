c ********** Capture Efficiency *******
c *** This subroutine calculates how the capture (aka retention)
c *** efficiency probability for each predator varies with prey size

       subroutine CaptureEfficiency(zmax, pnum,fnum,
     &         minind,maxind, zrmx,pr,capexp, 
     &         capture)


         integer zmax
         integer pnum
         integer fnum
         real*8 minind(fnum)
         real*8 maxind(fnum)
         real*8 zrmx(zmax,fnum)
         real*8 pr(pnum)
         real*8 capexp(fnum)
         real*8 capture(zmax,pnum,fnum)

         integer i,j,k
         integer ind1, ind2

         do 100 k = 1,fnum
             ind1 = minind(k)
             ind2 = maxind(k)
            do 200 i = ind1,ind2             
              do 300 j = 1,pnum
        
          capture(i,j,k) = ( 1-pr(j)/(zrmx(i,k)+pr(j)) )**capexp(k)

300            continue
200          continue
100      continue
         

       end ! subroutine
