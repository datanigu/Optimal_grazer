******  Prey physiological parameters *********
c *** Use this code to calculate phytoplankton physiological parameters
c *** having a choice between lienar or unimodal growth
c *** This code is based on "pparam.f"

c *** Calculating prey physiological parameters
      subroutine pparamchoice(pnum,pr,pv, mu0,muexp,muE,mutmpcoeff,
     &       kBoltz,tmp,pi,
     &       ks0,ksexp,etatmp, mu, ks, Dp, pm,u,muflag)
    
      implicit none
      integer pnum
      real*8 pr(pnum)
      real*8 pv(pnum)
      real*8 mu0
      real*8 muexp
      real*8 muE
      real*8 kBoltz
      real*8 tmp
      double precision pi
      real*8 mu(pnum)
      real*8 pm(pnum)
      real*8 u(pnum)
      real*8 ks0
      real*8 ksexp
      real*8 ks(pnum)
      real*8 etatmp
      real*8 Dp(pnum)
      real*8 mutmpcoeff
      integer muflag
      integer i

      do 45 i = 1,pnum

      ! Growth !
         if (muflag .eq. 1) then ! linear growth rate
          mu(i) =  (( mutmpcoeff*( (2.0*pr(i)) *10.0**6)**muexp )) 
     &      *exp(-muE/kBoltz/(273.15+tmp)) ! temperature-dependent growth (if include this last exponent and also a coefficient, but only matters if tmp does NOT equal 15 degrees C

c           mu = (( mu0*( (2.0*pr) *10.0**6)**muexp )) !non-temperature-dependent growth



          else if (muflag .eq. 2) then ! modal growth rate

              if (pv(i)*10.0**(18.0) .le. 40.0) then 
              mu(i) = !mutmpcoeff*(
     &      (10.0**(-0.43))*((pv(i)*10.0**18)**(0.19))/24.0/3600.0!)
!     &      *exp(-muE/kBoltz/(273.15+tmp)) ! temperature-dependent growth if include tecoeff and exponent

             else if (pv(i)*10.0**(18.0) .gt. 40.0) then
             mu(i) =! mutmpcoeff*(
     &        (10.0**(0.22))*((pv(i)*10.0**18)**(-0.15))/24.0/3600.0!)
!     &      *exp(-muE/kBoltz/(273.15+tmp)) ! temperature-dependent growth if include coeff and exponent
                end if!piecewise equation for unimodal growth

          end if !muflag



      ! Nutrient uptake !
      ks(i) = ks0*( (2.0*pr(i)) *10.0**6)**ksexp

      ! Viscosity !
      Dp(i) = ((tmp+273.15)*(1.38e-20))/(6.0*pi*etatmp*pr(i))

      ! Mortality !      
      pm(i) = pr(i)*0.0 + 0.005/24.0/3600.0

c     Swimming 
      u(i) = 0.0 + pr(i)*0.0

c      print*,'pm=',pm
c      print*,'u=',u
c      print*,'mu=',mu
c      print*,'ks=',ks
c      print*,'etatmp', etatmp
c      print*,'pi',pi
c      print*,'tmp',tmp
c      print*,'Dp=',Dp     
c
45      continue
      return
      end



