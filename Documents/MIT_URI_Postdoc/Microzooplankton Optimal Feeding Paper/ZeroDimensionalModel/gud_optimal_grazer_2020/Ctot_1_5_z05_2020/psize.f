c       Prey size class values 

      subroutine psize(pnum, psizewidth, psizemin, pi, pr,pv)
       implicit none
        double precision pi
        integer pnum

        real*8 psizewidth !width of size class in log space, usually 2 or 10
        real*8 psizemin !min size class, m
        real*8 pr(pnum) !phyto radius values, m
        real*8 pv(pnum) !phyto volume, m^3
        integer ii


c *** Calculating prey size classes
       do 200 ii = 1,pnum
          pr(ii) = psizemin*(psizewidth**(ii-1))
          pv(ii) = 4.0/3.0*pi*(pr(ii)**3)

200    continue

c       print*,'pv=',pv
 

       return
       end


