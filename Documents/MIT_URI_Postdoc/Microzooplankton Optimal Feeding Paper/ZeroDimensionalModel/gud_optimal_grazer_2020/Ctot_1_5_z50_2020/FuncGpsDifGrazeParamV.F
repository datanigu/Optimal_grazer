c ******* Functional Groups of Grazers *******

c **** This is the main code to examine change in prey and grazer biomass
c *** with time, including references to subprograms that determine 
c *** all the necessary parameters
c ***  This main code includes parameterizations for specific grazer types,
c ******namely dinoflagellates, ciliates, and a generic microzooplankton
c *** This code deals with the potential issue of really, really small biomass values
c *** as well as only records limited numbers of time steps
c *** This code is based on the code "FuncGpsExperimenting.f"
c *** Additionally, this code looks at the influence of different ways
c *** of describing the encounter kernel as well as variable capture probability
c *** Also, this code deals with changes in swimming speed depending on grazing intake
c ** This code also includes a change in handling time so that it is more realistic and may or may not have threshold
c *** This code is based on "FuncGpsDifParamIII.F" but now includes the correct grazing function (capture probability was in wrong place in earlier grazing versions) as well as a constant total carbon

c *** Subroutines included in this code are:
c ***** psize.f
c ***** pCcount.f
c ***** pparamchoice.f<--inludes monotonic or unimodal phyto growth
c ***** zsize.f
c ***** zCcountTypes.f
c ***** zparamtypes.f<--includes parameterization for dif. grazers, basal resp. corrected
c ***** handlingIII.f<--includes more reasonable values, with or without threshold
c ***** EncounterKernelTypes.f<--includes different kernels, including empirical ones
c ***** CaptureEfficiency.f
c ***** GrazingRateChoiceIV.f<--flexibility in how grazing is determined, capture probability in correct place
c ***** pzchangeII.f
c ***** swimspeed.f




c ** directory for this code is: cd /Users/dtaniguchi/Documents/MIT_URI_Postdoc/Grazer_Community_from_First_Principles/Functional_Groups_Fortran_Code/FuncGpsTypes

c *** This is what should type into terminal (since I still cannot figure out makefiles):
c gfortran pparamchoice.f psize.f pccount.f zCcountTypes.f zsize.f zparamtypes.f EncounterKernelTypes.f handlingIII.f GrazingRateChoiceIV.f pzchangeII.f ClearanceRate.f captureEfficiency.f swimspeed.f funcgpsdifgrazeparamV.f -o Ctot_1_5_z50_2018.exe

c *** Variables are declared in the .h file "FuncGpsParamTypes.h" 


      program FuncGpsDifGrazeParamIV

c *********************************************************************
c *** Declaring some local variables
c *********************************************************************

      implicit none
        double precision pi
        parameter (pi = 3.14159265359)
        integer tnum !number of time steps
        parameter(tnum = 1e7)
        real*8 dt !time step, s
        integer tcount !number to determine if should keep time step or not
        integer savestep ! how often want time changing variables to record

c ** Prey related
        integer pnum!number of prey size classes
        parameter(pnum=10)

c ** Grazer related
      integer fnum ! number of grazer functional groups
      parameter(fnum = 1)
      integer zmax !absolute max number of size classes, just a placeholder   
      parameter(zmax = 9)
      integer zfnum
      parameter(zfnum = zmax*fnum)! To initialize values to some number
 
c ** Other
      integer pzfnum
      parameter(pzfnum = pnum*zmax*fnum)
      integer pznum
      parameter(pznum = pnum*zmax)
      integer zftnum
      parameter(zftnum = zmax*tnum*fnum)
      integer ptnum
      parameter(ptnum = pnum*tnum)
      integer writeflag !set to 1 if want to write binary variables at end
      integer grazeflag ! determines how grazing is calculated
      integer encounterflag ! determines what sizes used for encounter kernel
      integer captureflag ! determines whether or not there is a capture 
                          ! (aka retention or catchability) efficiency probability
      integer mortalityflag ! determines if have linear or quadratic grazer mortality
      integer swimregflag(fnum)! determines if swimming determined by regression (=1)
c                 or based on multiple of body lengths (=2)
      integer swimflag(fnum) ! determines if swim speed is constant (=1) or
                            ! changes with grazing intake (=1)
      integer muflag ! determines if phytoplankton growth is unimodal or linear
      integer hflag ! determines how handling time is calculated
      integer debugflag !determines if will have print statements or not
      integer choosezflag ! set to 1 if want to define grazer size classes
 

c ** Indices for do loops
      integer t
      integer i
      integer kk
      integer ij
      integer ii
      integer pp
      integer ppz
      integer zz
      integer ff
      integer iiii
      integer ssi, ssii, mi,ma
      integer ind1, ind2,ind3
      integer indI, indII, indIII
      integer indo, indt
      integer indiv,indv
      integer q, st
      integer qq, qqq
      integer indq1, indq2
      
c ********* Referencing .h file for parameter declarations *******
      include "FuncGpsParamTypes.h"!parameters for main program and subprograms

c *** Flags
      writeflag = 1  ! set to 1 if want to write variables to binary files
      grazeflag = 1 !set flag to :
c                 1 if want Michaelis-Menten derivation to determine preferred grazing
c                 2 if want to use clearance rate to determine peferred grazing
c                 3 if want to impose a preferred prey size class (e.g., no switching)
      encounterflag = 5 !set to
c                 1 if want encounter kernel to just be based on prey
c                 2 if want encounter kernel to dynamically change based
c                    on size of prey relative to predator
c                 3 if want encounter kernel to be based on predator and prey size
c                       using Kiorboe and Titelman (1998, pg. 1617)
c                 4 if want encounter kernel based on empirical relationship
c                     To specify which relationship, use functype--
c                    functype = 1 is for dinoflagellates
c                    functype = 2 is for ciliates
c                    functype = 3 is for generic microzooplankton
c                 5 if want encounter kernel just based on predator and prey size
c                    
 
      captureflag = 2 ! set to:
c                 1 if want to have 100% capture efficiency for all
c                 2 if want to have capture efficiency decrease with incr. prey size 
      mortalityflag = 1! set to:
                        ! 1 if want linear grazer mortality
                        ! 2 if want quadratic grazer mortality
      swimregflag(1) = 2! set to 1 if want based on regression
c                               2 if want based on multiple of bodylengths
c      swimregflag(2) = 2


      swimflag(1) =1 ! set to:
                    ! 1 if want constant swim speed
                    ! 2 if want swim speed to change with grazing intake
c      swimflag(2) = 2


      muflag = 1! set to:
c                  1 if want linear change in phytoplankton growth with size
c                  2 if want unimodal change in phytoplankton growth with size
    
       hflag = 3 !set to
c            1 for handling based on gmax
c            2 for handling based on estimate of more reasonable times
c            3 for reasonable handling times but no threshold
c 
       debugflag = 0! set to
c            0 if don't want print statements
c            1 if want print statements for various variables in the code

       choosezflag = 0! set to 1 if want to define grazer size classes

c ********* Time **************
           dt = 100.0
           savestep = 10
           tcount = 1
c          print*,'dt',dt
c          numtemp = zmax*tnum*fnum


c ******* Defining parameters for environmental properties ***********
c (I don't think any of these need to be in matrix form because
c they get used in other variables which will become matrices)
      tmp = 15.0 !degrees C
      sal = 0.00035 !kg/kg = 35 psu
      etaW = 4.2844*10**(-5.0)+(0.157*(tmp+64.993)**2.0-91.296)**(-1.0)
      etatmp = etaW* 
     & ( 1+ (1.541+(1.998*10**(-2.0))*tmp-(9.52*10**(-5.0))*tmp**2.0)* 
     &     (sal) + 
     &     ((7.974-7.561*10**(-2.0))*tmp+
     &         (4.724*10**(-4.0))*tmp**2.0)*(sal**2.0))*1000!g/m/s
      kBoltz = 8.617*10.0**(-5)! eV/K

      ! Total carbon in system
      Ctot = 1.5!0.7 !gC/m^3
      !Cdis =  0.2!1.4! gC/m^3


c ******* Some grazer parameters ***********
c     functype is the type of functional group for the parameterization
c         1 = dinoflagellate
c         2 = ciliate
c         3 = generic microzooplankton
          functype(1) = 3
c          functype(2) = 2

c Min and max sizes for each functional group
      sizemin(1) =  25.0e-6 ! 25.0e-6! 250.0e-6! min size for first functional group
c      sizemin(2) = 10.0e-6! min size class for second functional group

      sizemax(1) = 100.0e-6! 100.0e-6! 1000.0e-6! approximate max size
c      sizemax(2) = 50.0e-6


c Base width of size class in log space (usually 2 or 10)
      zsizewidth = 10.0

c Carbon-associated values
      zCcoeff(1) = 0.3!0.1 !0.01
c      zCcoeff(2) = 0.01 
 
      zCexp(1) = 0.1
c      zCexp(2) = 0.1

c For encounter kernel (only important if encounterflag = 2)
      alpha = 0.0

c Physiological parameters (mortality, respiration, swimming, grazing)
      zgamma = 0.7
      zm0(1) = 0.05/24.0/3600.0   
c      zm0(2) = 0.05/24./3600
c      zm0(3) = 0.05/24./3600
      RE = 0.63!eV
      efficiency = 0.01 !efficiency of ciliary and flagellar motion

      gE = 0.67!eV, from Chen et al., 2012
      gmax0 = 33.96
      gmaxexp = -0.66
      gtmpcoeff = gmax0*exp(gE/(273.15+tmp)/kBoltz)

       capexp(1) = 3.0
c       capexp(2) = 3.0

       bodylengths(1) = 20.0!only important if swimregflag =2
c       bodylengths(2) = 10.0

       chi (1) = 0.1 !min swim speed is this fractional value of max swim speed
c       chi(2) = 0.1

       prefpreyclass(1) = 4 !integer that determines what prey to consume
c           when it doesn't dynamically change (e.g., grazeflag = 3)



c ************ Initializing some things to be zero ************
        data zchange/zftnum*0.0/
        data pchange/ptnum*0.0/
        data Cdischange/tnum*0.0/
        data gcountmx/pzfnum*0.0/
        data gtimechange/zftnum*0.0/
        data clearance/pzfnum*0.0/
        data beta/pzfnum*0.0/
        data v/zfnum*0.0/! Initializing all swimming speeds as zero
        Cdischangetemp = 0.0


c *********************************************************************
c *********** Prey-related parameters and variables ************ 
c *********************************************************************

c *** Values for prey size parameters
       psizewidth = 2.0
       psizemin = 0.5e-6

c *** Calculating prey size classes
       call psize(pnum, psizewidth, psizemin,pi, pr, pv)
     
c *** Calculating prey carbon content and initial number of and 
c *** carbon concentration in each size class 
       pCcoeff = 0.1
       pCexp = 0.1
       call pCcount(pnum,pr,pv,pCcoeff,pCexp,pCsize,pCcell,pcount)

       pchange(:,1) = pCsize
       pchangetemp1 = pCsize !initializing temporary variable


c *** Calculating physiolgoical parameters
      ! Growth !
      mu0 = 1.36/24.0/3600.0!1/s
      muexp = -0.16
      muE = 0.36 ! from Chen et al., 2012
c     mutmpcoeff  = 1.5804e6 ! average value from temperature correction in Taniguchi et al., 2014
      mutmpcoeff = mu0*exp(muE/(tmp+273.15)/kBoltz)

       if (debugflag .eq. 1) then
      print*,'pr*10^6',pr*(10.0**6)
       endif

      ! Nutrient uptake !
      ks0 = 0.33*106.0/16.0*12.0/1000.0!gC/m^3
      ksexp=0.48
          
c      Use pparam if do NOT want a choice in phytoplankton growth shape
c      call pparam(pnum,pr, mu0,muexp,muE,mutmpcoeff,kBoltz,tmp,pi, 
c     &       ks0,ksexp,etatmp, mu, ks, Dp, pm,u)
          

c **     Use pparamchoice if want to choose shape of size-depend. of phyto. growth
      call pparamchoice(pnum,pr,pv,mu0,muexp,muE,mutmpcoeff,
     &    kBoltz,tmp,pi, 
     &       ks0,ksexp,etatmp, mu, ks, Dp, pm,u,muflag)


c       print*,'pr',pr
c       print*,'pv',pv
c       print*,'mu=',mu
c       print*,'ks = ',ks

c *********************************************************************
c *********** Grazer-related parameters and variables **************
c *********************************************************************
c *** Calculating grazer size classes
      call zsize(zmax,fnum,pi,sizemax,sizemin,zsizewidth,
     &     zrmx,znumf,zvmx,minind,maxind)

c         print*,'zrmx*10^6',zrmx*10.0**6.0
c

c         print*,'minind outside subroutine',minind
c         print*,'maxind outside subroutine',maxind

c      if (choosezflag .eq. 1) then
c         zrmx(1,1) = 5.0e-6
c         zrmx(2,1) = 25.0e-6
c         zrmx(3,1) = 50.0e-6
c      endif
c
         print*,'zrmx*10^6',zrmx*10.0**6.0
c          print*,'zvmx',zvmx
c         print*,'zrall',zrall
       
      call  zCcountTypes(fnum,zmax,zrmx,zvmx,zCcoeff,zCexp,
     &     functype,
     &     zCcell,zCsize,zcount)

c       print*,'zCsize',zCsize
        do 220 kk = 1,fnum
           zchange(:,1,kk) = zCsize(:,kk)
           zchangetemp1(:,kk) = zCsize(:,kk)!temporary variable
220     continue

c         print*,'zchange(:,1,1)',zchange(:,1,1)
c         print*,'zchange(:,1,2)',zchange(:,1,2)


c       print*,'gE',gE
c       print*,'gmax0',gmax0
c       print*,'exp(gE/(273+tmp)/kBoltz)',exp(gE/(273.15+tmp)/kBoltz)
c      print*,'zrmx',zrmx

c *** Calculating pysiological parameters
        call zparamtypes(fnum,zmax,znumf,zrmx,zvmx,pi,
     &   minind, maxind,tmp,etatmp,
     &   gtmpcoeff, gmaxexp,gE,kBoltz,
     &   bodylengths,functype,
     &   efficiency,zm0,RE,zCcell,
     &   v,Dz,gmax,zm,R0,Rv,R,swimregflag)

        do 132 ssi=1,fnum
            mi = minind(ssi)
            ma = maxind(ssi)
          do 133 ssii = mi,ma
          
              vmax(ssii,ssi) = v(ssii,ssi)
              vmin(ssii,ssi) = chi(ssi)*v(ssii,ssi)
133       continue
132    continue
       
          if (debugflag .eq. 1) then
c        print*,'vmin',vmin

c         print*,'vmax',vmax        


       
c       print*,'functype',functype
c       print*,'swimregflag',swimregflag
c       print*, 'v=',v
c       print*,'zvmx',zvmx
c       print*, 'Dz=',Dz
c       print*, 'zm=',zm
c       print*, 'R=',R
c       print*, 'zgamma',zgamma
c       print*, 'zm0',zm0
c       print*, 'R0 right after zparamtypes',R0
c       print*, 'efficiency',efficiency
c       print*, 'gmax',gmax
c       print*, 'R', R
c
c        print*,'gmax',gmax
c        print*,'gmaxexp',gmaxexp
c        print*,'gmax(:,1)   ',gmax(:,1)
c        print*,'gmax(:,2)   ',gmax(:,2)
          endif

c *** Calculating dissolved carbon
      Cdis = Ctot - sum(zCsize) - sum(pCsize)
       if (Cdis .lt. 0.0) then
         print*,'Negative dissolved carbon'
       endif

c      print*,'Ctot',Ctot
      Cdischange(1) = Cdis
      Cdischangetemp = Cdis !initializing temporary variable
 
c      print*,'Initial Cdishchange(1)',Cdischange(1)
c      print*,'Initial Cdischangetemp',Cdischangetemp

c *** Calculating handling time
        
       call handlingIII (zmax, pnum,fnum,znumf,
     &   gmax, pCcell, zCcell,zrmx,pr,
     &    minind,maxind, hflag,
     &     h)

c       call handling (zmax, pnum,fnum,znumf,
c     &   gmax, pCcell, zCcell,
c     &    minind,maxind, 
c     &     h)
c

       if (debugflag .eq. 1) then
       print*,'h(:,2,1)',h(:,2,1)
       endif

c *** Calculating encounter kernel
       call EncounterKernelTypes (zmax, pnum,fnum,znumf,
     &     Dp, Dz, pr, zrmx,zvmx, pi,minind,v,alpha,functype,
     &     beta,encounterflag)

c       print*,'beta(:,5,1)',beta(:,5,1)
c       print*,'beta(:,1,2)',beta(:,1,2)
c
c      print*,'zrmx*10^6',zrmx*10.0**6
c      print*,'log beta(:,1,1)',log(beta(:,1,1))
c      print*,'log beta(:,2,1)',log(beta(:,2,1))
c      print*,'log beta(:,3,1)',log(beta(:,3,1))
c      print*,'log beta(:,4,1)',log(beta(:,4,1))
c      print*,'log beta(:,5,1)',log(beta(:,5,1))
c      print*,' '
c      print*,'log beta(:,1,2)',log(beta(:,1,2))
c      print*,'log beta(:,2,2)',log(beta(:,2,2))
c      print*,'log beta(:,3,2)',log(beta(:,3,2))
c      print*,'log beta(:,4,2)',log(beta(:,4,2))
c      print*,'log beta(:,5,2)',log(beta(:,5,2))
c

c *** Calculating clearance rate

        call ClearanceRate(zmax,fnum,pnum,minind, maxind,
     &    gmax, zCcell,beta,h,clearance)

c        print*,'clearance(:,:,1)',clearance(:,:,1)
c        print*,'clearance(:,:,2)',clearance(:,:,2)

c *** Calculating capture efficiency (aka retention, catchability)
      if (captureflag .eq. 1) then ! able to retain everything with 100% efficiency
       data capture/pzfnum*1.0/ 

       else !Calculate change in capture efficiency with size of prey
       data capture/pzfnum*0.0/

      call CaptureEfficiency(zmax, pnum,fnum,
     &         minind,maxind, zrmx,pr,capexp, 
     &         capture)
        endif


c      print*,'capture(:,1,1)',capture(:,1,1)
c      print*,'capture(:,2,1)',capture(:,2,1)
c      print*,'capture(:,3,1)',capture(:,3,1)
c      print*,'capture(:,4,1)',capture(:,4,1)
c      print*,'capture(:,5,1)',capture(:,5,1)
c      print*,' '
c      print*,'capture(:,1,2)',capture(:,1,2)
c      print*,'capture(:,2,2)',capture(:,2,2)
c      print*,'capture(:,3,2)',capture(:,3,2)
c      print*,'capture(:,1,2)',capture(:,4,2)
c      print*,'capture(:,1,2)',capture(:,5,2)
c 


c***********************************************************************
c ********* Beginning loop to calculate change in biomass with time ****
c***********************************************************************

        do 3003 t = 1,tnum-1

c *** Initializing some things to zero

        do 101 pp = 1,pnum
         gprey(pp) = 0.0
          do 202 ppz = 1,zmax
           gpressure(ppz,pp) = 0.0
           ggrazerbio(ppz,pp) = 0.0
           gtemp(ppz,pp) = 0.0
202       continue
101     continue

        do 303 ff = 1,fnum
           do 404 zz = 1,zmax
             ggrazer(zz,ff) = 0.0
             gtimetemp(zz,ff) = 0.0       
404         continue
303     continue
       Ctemp = 0.0

c
c       data gprey/pnum*0.0/
c       data ggrazer/zfnum*0.0/
c       data gtimetemp/zfnum*0.0/
c       data gpressure/pznum*0.0/
c       data ggrazerbio/pznum*0.0/
c       data gtemp/pznum*0.0/

c           print*,'t',t
c           print*,' '
cc           print*,' '
c           print*,'gprey new time',gprey
c           print*,'Ctemp',Ctemp
c           print*,'gtimetemp',gtimetemp
c 
c *** Calculating which size class grazers consume
 
c      print*,'znumf(ii)',znumf(ii) 
c      print*,'ind1',ind1
c      print*,'ind2',ind2
c      print*,'ind3',ind3
c      print*,'ggrazer(ind1:ind2,1) in gr 1', ggrazer(ind1:ind2,1)
c       print*,'testing',testing
c        print*,'gprey before gr',gprey
c        print*,'zchange(ind1:ind2,t,ind3)',zchange(ind1:ind2,t,ind3)
c        print*,'pchangetemp1',pchangetemp1
c        print*,'pchangetemp2',pchangetemp2
c      print*,'zchangetemp1(ind1:ind2,ind3)',zchangetemp1(ind1:ind2,ind3)
c       print*,'zchangetemp1(:,1)',zchangetemp1(:,1)
c       print*,'zchangetemp1(:,2)',zchangetemp1(:,2)
c           print*,' '
c           print*,' '
c



cc        print*,'grazeflag before subroutine',grazeflag

c         print*,'minind outside subroutine',minind
c         print*,'maxind outside subroutine',maxind



        call  GrazingRateChoiceIV (pnum,zmax,fnum,! Choice in grazing preference
     &   minind,maxind,
     &   pCcell, zCcell, 
     &   pchangetemp1,
     &   zchangetemp1,
     &   beta,h,gmax,clearance,
     &   capture,grazeflag, prefpreyclass,  
     &   gtimetemp,gtemp,
     &   gprey,ggrazer)


c      print*,'ggrazer(ind1:ind2,1) in gr 2', ggrazer(ind1:ind2,1)
c       print*,'testing',testing

c       ggrazer(ind1:ind2,ind3) = testing(1:znumf(ii),ii)
      
c       gcountmx(ind1:ind2,:,ii) = gcountmx(ind1:ind2,:,ii) + gcount

c       print*, 'gtimetemp',gtimetemp
c       print*,'gtimechange(:,t,ii)',gtimechange(:,t,ii)

c      print*,'gprey after grazing rate choice',gprey   
c      print*,'ggrazer(:,1) after grazing rate choice',ggrazer(:,1)
c      print*,'ggrazer(:,2) after grazing rate choice',ggrazer(:,2)
c
c *** Print check for ggrazer
c      print*,'ggrazer(ind1:ind2,1) out gr', ggrazer(ind1:ind2,1)
c      print*,'ggrazer(ind1:ind2,2) out gr', ggrazer(ind1:ind2,2)
c      print*,'testing outside',testing



c ** Finding biomass change
c           print*,'iii',iii 
c           print*,'testing iii init',testing
c           print*,'ggrazer',ggrazer
c        print*,'ggrazer(:,iii)in iii loop',ggrazer(indI:indII,indIII)
c        print*,'zchange(ind1:ind2,t,ind3)',zchange(indI:indII,t,ind3)


c**     pzchange does not allow potential for quadratic mortality closure 
c          call pzchange(pnum,zmax,fnum,minind,maxind, 
c     &     dt,Ctot,
c     &     mu, Cdischangetemp, ks, pm, gprey,
c     &     zm, R,
c     &     zgamma,ggrazer,
c     &     pchangetemp1, pchangetemp2,
c     &    zchangetemp1,
c     &   zchangetemp2)

c**      pzchangeII allows for potential quadradic mortality for predator
          call pzchangeII(pnum,zmax,fnum,minind,maxind, 
     &     dt,Ctot,
     &     mu, Cdischangetemp, ks, pm, gprey,
     &     zm, R,
     &     zgamma,ggrazer,
     &     mortalityflag,
     &     pchangetemp1, pchangetemp2,
     &    zchangetemp1,
     &   zchangetemp2)


c       print*,'testing iii fin',testing
c        print*,'pchangetemp1',pchangetemp1
c        print*,'pchangetemp2',pchangetemp2
c      print*,'zchangetemp1(ind1:ind2,ind3)',zchangetemp1(ind1:ind2,ind3)


c *** Determining swim speed based on grazing intake ******************

      do 124 qq = 1,fnum
       if (swimflag(qq) .eq. 2) then 

          indq1 = minind(qq)
          indq2 = maxind(qq)
         do 125 qqq =indq1,indq2
            g(qqq,qq) = ggrazer(qqq,qq)/zchangetemp1(qqq,qq)
125       continue
     

      call swimspeed(fnum,zmax,maxind,minind,
     &     vmax,vmin, R0, Rv, g,v)


c       print*,'time step, t=',t
c       print*,'v from swimspeed subroutine = ',v

c *** Calculating new encounter kernel, clearance rate, and respiration based on new swim speed
            call EncounterKernelTypes (zmax, pnum,fnum,znumf,
     &     Dp, Dz, pr, zrmx,zvmx, pi,minind,v,alpha,functype,
     &     beta,encounterflag)

        call ClearanceRate(zmax,fnum,pnum,minind, maxind,
     &    gmax, zCcell,beta,h,clearance)

           do 122 q=1,fnum
               indiv = minind(q)
               indv = maxind(q)
            do 123 st= indiv,indv
      Rv(st,q) = ((3.0*pi*zrmx(st,q)*2.0*(v(st,q)**2.0))*etatmp/1000.0
     &        *(10.0**(-3.0))/20.2/(tmp+273.15) 
     &        / 0.082*12.0/zCcell(st,q)/efficiency)

123         continue
122        continue

      R = Rv + R0     

        endif !swimflag
124    continue  


c *** Updating the dissolve nutrient value,
c ** grazer biomass, and which size class grazed 
       do 2002 ii=1,fnum
          ind1 = minind(ii)! defining some indices
          ind2 =  maxind(ii)
       gtimechange(ind1:ind2,t,ii) = gtimetemp(ind1:ind2,ii)!should change each time
       Ctemp = sum(zchangetemp2(ind1:ind2,ii)) + Ctemp
       zchangetemp1(ind1:ind2,ii) = zchangetemp2(ind1:ind2,ii)
2002    continue

        pchangetemp1 = pchangetemp2
        Cdischangetemp = Ctot - sum(pchangetemp2) - Ctemp



c ************ Only saving some time steps *************
                 if (mod(t, savestep) .eq. 0) then
c         print*,'t',t
c         print*,'mod(t, savestep)',mod(tnum, savestep) 

         tcount = tcount + 1
         pchange(:,tcount) = pchangetemp2
c         Cdischange(tcount) = Ctot - sum(pchangetemp2) - Ctemp
          Cdischange(tcount) = Cdischangetemp

c          print*,'Cdischange(tcount)',Cdischange(tcount)
c         print*,'tcount',tcount

          do 4004 iiii = 1,fnum
            indo = minind(iiii)
            indt = maxind(iiii)
           zchange(indo:indt,tcount,iiii) = 
     &  zchangetemp2(indo:indt,iiii)
4004     continue

                   end if ! for modular arthimetic
c******************************************************************

            

c *** Udating other time changing variables
c       Cdischange(t+1) = Ctot - sum(pchange(:,t+1))-Ctemp! Only use for all timesteps

c           print*,' '
c           print*,'t',t
c           print*,'Cdischangetemp',Cdischangetemp

3003     continue ! end time stepping loop

c        pchange = pchange(:,1:tcount)
c        zchange = zchange(:,1:tcount,:)
c        Cdischange = Cdischange(1:tcount)
c
c
c ** Print check **
        if (debugflag .eq. 1) then
c       print*,'zrmx*10^6=',zrmx*10**6
c       print*,'pr*10^6=',pr*10**6


        print*,'pchange(:,1)',pchange(:,1)
        print*,'pchange(:,2)',pchange(:,2)
        print*,'pchange(:,3)',pchange(:,3)

        print*,'pchange(:,tcount-1)',pchange(:,tcount-1)
        print*,'pchange(:,tcount)',pchange(:,tcount)
        print*,'log(pchange(:,tcount-1))',log(pchange(:,tcount-1))
        print*,'log(pchange(:,tcount))',log(pchange(:,tcount))


        print*,' '
        print*,'zchange(:,1,1)',zchange(:,1,1)
        print*,'zchange(:,2,1)',zchange(:,2,1)
         print*,'zchange(:,3,1)',zchange(:,3,1)

        print*,' '

       print*,'zchange(:,tcount-1,1)',zchange(:,tcount-1,1)
        print*,'zchange(:,tcount,1)',zchange(:,tcount,1)
        print*,' '

        print*,'log(zchange(:,1,1))',log(zchange(:,1,1))
  
        print*,' '

 
        print*,'log(zchange(:,tcount-1,1))',log(zchange(:,tcount-1,1))
        print*,'log(zchange(:,tcount,1))',log(zchange(:,tcount,1))


        print*,' '
        print*,' '
        print*,' '

c        print*,'zchange(:,1,2)',zchange(:,1,2)
c        print*,'zchange(:,2,2)',zchange(:,2,2)
c        print*,'zchange(:,3,2)',zchange(:,3,2)
c
c        print*,'zchange(:,tcount-1,2)',zchange(:,tcount-1,2)
c        print*,'zchange(:,tcount,2)',zchange(:,tcount,2)
c        print*,' '
c        print*,'log(zchange(:,1,2))',log(zchange(:,1,2))
c 
c        print*,' '
c        print*,'log(zchange(:,tcount-1,2))',log(zchange(:,tcount-1,2))
c        print*,'log(zchange(:,tcount,2))',log(zchange(:,tcount,2))
c


        print*,' '
        print*,'Cdischange(1)',Cdischange(1) 
        print*,'Cdischange(2)',Cdischange(2) 

        print*,'Cdischange(tcount-1)',Cdischange(tcount-1) 
        print*,'Cdischange(tcount)',Cdischange(tcount)


cc        print*,' '
cc        print*,'gtimechange(:,1,1)',gtimechange(:,1,1)
cc        print*,'gtimechange(:,tcount-1,1)',gtimechange(:,tcount-1,1)
cc        print*,'gtimechange(:,tcount,1)',gtimechange(:,tcount,1)
cc
cc        print*,' '
cc        print*,'gtimechange(:,1,2)',gtimechange(:,1,2)
cc        print*,'gtimechange(:,tcount-1,2)',gtimechange(:,tcount-1,2)
cc        print*,'gtimechange(:,tcount,2)',gtimechange(:,tcount,2)
cc

        print*,'R0',R0
        print*,'R',R
        print*,'beta',beta
        print*,'capture',capture

        print*,'dt',dt
        print*,'tnum',tnum
        print*,'tcount',tcount
c
c
c
c



c         print*,'dt',dt
c         print*,'zrmx',zrmx
c        print*,'Ctot',Ctot
c         print*,'gcountmx(1,:,1)',gcountmx(1,:,1)
c        print*,'gcountmx(2,:,1)',gcountmx(2,:,1)
c        print*,'gcountmx(3,:,1)',gcountmx(3,:,1)
c        print*,'gcountmx(4,:,1)',gcountmx(4,:,1)
c        print*,'gcountmx(5,:,1)',gcountmx(5,:,1)
c
c         print*,'gcountmx(:,1,1)',gcountmx(:,1,1)
c         print*,'sum(gcountmx(:,1,1)',sum(gcountmx(:,1,1))
c         print*,'gtimechange(:,1,1)',gtimechange(:,1,1)
c         print*,'gtimechange(:,1,2)',gtimechange(:,1,2)
          

          endif ! debugflag 
c
c ********************************************************************
c ************** Writing variables of interest to binary form ********
        if (writeflag .eq. 1) then

c ** Biomass changes with time variables ***
c Prey
        reclen = pnum*(tcount)*8
        open(17,file='pchange.bin',status='replace',access='direct',
     &  recl = reclen)
        write(17,rec=1)pchange(:,1:tcount)
        close(17)

c Grazers
        reclen = zmax*fnum*tcount*8
        open(17,file='zchange.bin',status='replace',access='direct',
     &  recl = reclen)
        write(17,rec=1)zchange(:,1:tcount,:)
        close(17)

c Dissolved and total carbon
        reclen = tcount*8
        open(17,file='Cdischange.bin',status='replace',
     &  access='direct', recl = reclen)
        write(17,rec=1)Cdischange(1:tcount)
        close(17)

        reclen = 8
        open(17,file='Ctot.bin',status='replace',
     &  access='direct', recl = reclen)
        write(17,rec=1)Ctot
        close(17)

c ** Grazing count
       reclen = zmax*pnum*fnum*8
        open(17,file='gcountmx.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)gcountmx
         close(17)

      reclen = zmax*tnum*fnum*8!numtemp
        open(17,file='gtimechange.bin',status='replace',access='direct',
     & recl = reclen)
        write(17,rec=1)gtimechange
        close(17)


c ** Time variables
      reclen = 8
        open(17,file='dt.bin',status='replace',access='direct', 
     &  recl = reclen)
        write(17,rec=1)dt
        close(17)

        reclen = 8
        open(17,file='tnum.bin',status='replace',access='direct', 
     &  recl = reclen)
        write(17,rec=1)tnum
        close(17)

        reclen = 8
        open(17,file='tcount.bin',status='replace',access='direct', 
     &  recl = reclen)
        write(17,rec=1)tcount
        close(17)

        reclen = 8
        open(17,file='savestep.bin',status='replace',access='direct', 
     &  recl = reclen)
        write(17,rec=1)savestep
        close(17)



c ** Prey size values
       reclen = 8
        open(17,file='pnum.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)pnum
         close(17)

       reclen = pnum*8
        open(17,file='pr.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)pr
         close(17)

       reclen = pnum*8
        open(17,file='pCcell.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)pCcell
         close(17)

       reclen = pnum*8
        open(17,file='mu.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)mu
         close(17)


c ** Grazer size values
       reclen = 8
        open(17,file='zmax.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)zmax
         close(17)

       reclen = fnum*8
        open(17,file='znumf.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)znumf
 

       reclen = fnum*8
        open(17,file='minind.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)minind
         close(17)

       reclen = fnum*8
        open(17,file='maxind.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)maxind
         close(17)


       reclen = 8
        open(17,file='fnum.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)fnum
         close(17)

       reclen = zmax*fnum*8
        open(17,file='zrmx.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)zrmx
         close(17)

       reclen = zfnum*8
        open(17,file='zCcell.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)zCcell(:,1)
         close(17)

      reclen = 8
        open(17,file='zsizewidth.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)zsizewidth
         close(17)

C ** Swim speeds
       reclen = zfnum*8  
        open(17,file='v.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)v
         close(17)

c ** Flags
       reclen = 8
        open(17,file='grazeflag.bin',status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)grazeflag
         close(17)

       reclen = 8
        open(17,file='encounterflag.bin',
     &  status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)encounterflag
         close(17)

       reclen = 8
        open(17,file='captureflag.bin',
     &  status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)captureflag
         close(17)

       reclen = 8
        open(17,file='mortalityflag.bin',
     &  status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)mortalityflag
         close(17)

        reclen = 8
        open(17,file='muflag.bin',
     &  status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)muflag
         close(17)

     
       reclen = 8
        open(17,file='swimflag.bin',
     &  status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)swimflag
         close(17)


       reclen = 8
        open(17,file='swimregflag.bin',
     &  status='replace',access='direct',
     &  recl= reclen)
        write(17,rec=1)swimregflag
         close(17)

 




      endif!writeflag



      end ! Main program











