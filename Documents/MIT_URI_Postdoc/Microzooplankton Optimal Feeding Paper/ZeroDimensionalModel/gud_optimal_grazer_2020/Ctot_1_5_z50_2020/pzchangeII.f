c **********  Change in prey and grazer biomass with time *******

c *** Use this subroutine to calculate change in bioass wtih each time step


      subroutine pzchangeII(pnum,zmax,fnum,minind,maxind, 
     &     dt,Ctot,
     &     mu, Cdis, ks, pm,gprey,
     &     zm, R, 
     &    zgamma, ggrazer,
     &    mortalityflag,
     &     ptemp, ptemp2,
     &     ztemp, ztemp2)
c     &     Cdis2)
   

c *** Defining variables 
      implicit none
      integer zmax, pnum,fnum
      real*8 minind(fnum)
      real*8 maxind(fnum)
      real*8 dt !time step
      real*8 Ctot
      real*8 mu(pnum)
      real*8 Cdis
      real*8 ks(pnum)
      real*8 pm(pnum)
      real*8 gprey(pnum)
      real*8 zm(zmax,fnum)
      real*8 R(zmax,fnum)
      real*8 zgamma
      real*8 ggrazer(zmax,fnum)
      integer mortalityflag
      real*8 ptemp(pnum)
      real*8 ptemp2(pnum)
      real*8 ztemp(zmax,fnum)
      real*8 ztemp2(zmax,fnum)
c      real*8 Cdis2


      integer i, j,k


c      print*,'gprey',gprey
c      print*,'ggrazer',ggrazer

c ******* Print check
c      print*,' '
c      print*,' '
c      print*,' '
c
c      print*, 'Entering change in biomass with time subprogram'
c      print*,'ptemp',ptemp
c      print*, 'Cdis',Cdis
c      print*,'mu',mu
c      print*,'ks',ks
c      print*,'pm',pm
c      print*,'gprey',gprey
c      print*,' '
c      print*,' '
c      print*,'ztemp',ztemp
c      print*,'zm',zm
c      print*,'R',R
c      print*,'zgamma',zgamma
c      print*,'ggrazer',ggrazer
c      print*,'dt',dt
c
c      print*,' '
c      print*,' '


c ********** Governing equations of model **************** c

c ** Calculating biomass values at new time step
       do 130 i = 1,pnum
         ptemp(i) = max(ptemp(i), 10.0**(-20.0))! prevents ptemp from being negative

          ptemp2(i)  = ptemp(i)+ 
     &       (ptemp(i) * 
     &       (mu(i)*Cdis/(Cdis + ks(i)) - pm(i)) -
     &       gprey(i) )*dt

       ptemp2(i) = max(ptemp2(i),10.0**(-20.0))! prevents ptemp2 from being negative
130    continue


      if (mortalityflag .eq. 1) then ! linear grazer mortality
      do 140 k=1,fnum
          do 150 j = 1,zmax

          ztemp(j,k) = max(ztemp(j,k), 10.0**(-20.0))! prevents ztemp from being negative

c ********* change in grazer biomass, respiration constant (but based on cell size); based on Fenchel and Finlay (1983)
c          ztemp2(j,k) = ztemp(j,k) + 
c     &        (ztemp(j,k)*(-zm(j,k) - R(j,k)) + zgamma*ggrazer(j,k))*dt

c **** Respiratory loss scales with ingestion, based on Verity (1985)
          ztemp2(j,k) =  ztemp(j,k) + 
     &        (ztemp(j,k)*(-zm(j,k) - R(j,k))!+0.00000278 )  
     &        - (0.63*ggrazer(j,k))
     &        + zgamma*ggrazer(j,k))*dt

          ztemp2(j,k) = max(ztemp2(j,k), 10.0**(-20.0))! prevents ztemp2 from being negative

150       continue
140   continue

        else if (mortalityflag .eq. 2) then ! quadratic grazer mortality
         do 160 k=1,fnum
          do 170 j = 1,zmax

          ztemp(j,k) = max(ztemp(j,k), 10.0**(-20.0))! prevents ztemp from being negative

          ztemp2(j,k) = ztemp(j,k) + 
     &  (ztemp(j,k)*(-zm(j,k)*ztemp(j,k) - R(j,k))
     &   + zgamma*ggrazer(j,k))*dt

          ztemp2(j,k) = max(ztemp2(j,k), 10.0**(-20.0))! prevents ztemp2 from being negative

170       continue
160   continue

      end if

     

cc *** Updating the dissolve nutrient value
c      Cdis2 = Ctot - sum(ptemp2) - sum(ztemp2)

c*** Print check
c      print*,'ptemp2',ptemp2
c
c       print*,'ztemp2',ztemp2
c
c      print*,'Cdis2',Cdis2
c      print*,' '
c      print*,' '
c


cc           testp =  (mu(1,:)*Cdis/(Cdis + ks(1,:)) - pm(1,:)) -
cc     &       gprey 
cc           testz = ztemp(:,1)*(-zm(:,1) - R(:,1)) + zgamma*ggrazer
cc
ccc        print*,'testp',testp
ccc        print*, 'testz',testz
cc




c *** Print check
c             print*, 'gprey after biomass calc', gprey
c             print*, 'ggrazer after biomass calc',ggrazer
c             print*, 'gtemp after biomass calc',ztemp
c             print*, 'gtemp2 after biomass calc',ztemp2
c             print*, 'ptemp after biomass calc',ptemp
c             print*, 'ptemp2 after biomass calc',ptemp2
c             print*, 'Cdis after change',Cdis
c  



c
c      print*, 'Leaving change in biomass with time subprogram'

      end !pz_change subroutine



