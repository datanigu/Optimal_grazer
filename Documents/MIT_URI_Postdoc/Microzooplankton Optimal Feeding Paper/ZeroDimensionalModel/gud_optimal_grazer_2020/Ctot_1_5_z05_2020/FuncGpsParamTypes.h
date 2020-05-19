c ****** Parameters for grazer and prey groups and environmental parameters ********

c *************** Environmental and misc. parameters **********
      real*8 tmp! temperature, C
      real*8 etaW ! kinematic viscosity of pure water, kg/m/s
      real*8 etatmp !temperature-dep kinematic viscosity for respiration, g/m/s, 
                    ! from Sharqawy et al., 2010, equation 22
      real*8 sal !salinity, kg/kg
      real*8 kBoltz ! Boltzmann's contant = 8.617*10^-5 eV/K (where K = Kelvin)
      real*8 Cdis !dissolved carbon, gC/m^3
      real*8 Ctot !total carbon in system, gC/m^3


c *************** Prey parameters **************
c *** Parameters for size
        real*8 psizewidth !width of size class in log space, usually 2 or 10
        real*8 pr(pnum) !phyto radius values, m
        real*8 psizemin !min size class, m

c *** Parameters for carbon concentration
        real*8 pv(pnum)! cell volumne, m^3
        real*8 pCcoeff! coefficient for carbon per size class, g/m^3
        real*8 pCexp! dimensionless, exponent for carbon per size class
        real*8 pCsize(pnum)! carbon concentration in each size class, gC/m^3        
        real*8 pCcell(pnum)! carbon per cell, gC/cell
        real*8 pcount(pnum)! number of cells per vol per size class, cells/m^3

c *** Parameters for physiological rates
c ptemp = intermediate value for calculating change in prey biomass with time, gC/m^3
c ptemp2 = anther intermediate value, gC/m^3
c p_change = change in prey biomass with time, gC/m^3, 
c       where each row is a time point, and each column is a prey size class
      real*8 mu0 !prey growth coefficient, 1/s, taken from Taniguchi et al. 2014
      real*8 muexp ! exponent for prey growth, dimensionless
      real*8 mutmpcoeff ! coefficient that includes temperature-dep
      real*8 muE !activation energy for growth rate, eV (electron volts)
      real*8 mu(pnum) !prey growth (may or may not be temperature dependent), 1/s
      real*8 pm(pnum) ! prey general loss rate, 1/s
      real*8 ks0 !nutrient half saturation constant = ks0*diameter^ksexp, gC/m^3 
      real*8 ksexp !nutrient half saturation constant exponent
      real*8 ks(pnum) ! nutrient half saturation constant = ks0*diameter^ksexp, gC/m^3
      real*8 u(pnum) ! prey swimming speed, m/s
      real*8 Dp(pnum) ! diffusion coefficient

c *************** Grazer parameters ***********************
      integer functype(fnum)
c ****** Grazer size-related parameters
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
 
c *** Parameters for finding carbon-associated values **********
      real*8 zCcoeff(fnum)! coefficient for carbon per size class
      real*8 zCexp(fnum) ! exponent for carbon per size class
      real*8 zCsize(zmax,fnum)! carbon per size class, gC/m^3/size class
      real*8 zcount(zmax,fnum) !cells per size class, cells/zie class
      real*8 zCcell(zmax,fnum) !carbon per cell, gC/cell

c *** Parameters to find physiological parameters
      real*8 zgamma ! assimilation efficiency
      real*8 vmin(zmax,fnum) !min swimming speed for each func gp, m/s
      real*8 vmax(zmax,fnum) !max swim speed, m/s
      real*8 chi(fnum) ! min swim speed = chi*max swim speed, dimensionless
      real*8 v(zmax, fnum)
      real*8 vsizewidth !log base size width (if not linear), usually 2 or 10
      real*8 bodylengths(fnum)!mutliple of body lengths swim if not using regression
c      real*8 betav(zmax,fnum)
      real*8 alpha !scalar to determine transition between encounter kernel
c            based on just prey size or predator + prey size
      real*8 Dz(zmax,fnum) !diffusion coefficient for grazers
      real*8 gmax0 !max grazing coefficient includng temperature-dep
      real*8 gtmpcoeff ! max grazing coeff not including temp-dep
      real*8 gmaxexp !max grazing exponent 
      real*8 gE ! grazing activation energy, eV
      real*8 gmax(zmax,pnum,fnum) !max grazing rate, 1/s
      real*8 g(zmax,fnum) ! grazing rate not multiplied by grazer biomass, 1/s
      real*8 zm0(fnum) !basal mortality rate (i.e., not due to swimming), 1/s      
      real*8 zm(zmax,fnum)! 1/s, total grazing mortality, from swimming and basal 
      real*8 z0R(zmax,fnum) !O2/cell/hr, basal respiration for grazer
      real*8 R0(zmax,fnum) !1/s, specific basal respiration for grazer  
      real*8 RE ! eV, actiation energy for respiration
      real*8 efficiency !efficiency of ciliary/flagellar propulsion
      real*8 Rv(zmax,fnum) ! 1/s, specific respiration from swimming
      real*8 R(zmax,fnum) !1/s, total respiration (basal and swimming)  
      integer prefpreyclass(fnum) ! prey class to consume if 
c           want it to not dynamically change


c *** Parameters calculated from subroutines
      real*8 h(zmax,pnum,fnum) !handing time, s
      real*8 beta(zmax,pnum,fnum)
      real*8 gpressure(zmax,pnum)
c      real*8 gcount(zmax,pnum)
c      real*8 gtimetemp(zmax)
      real*8 gtimetemp(zmax,fnum)
      real*8 gtemp(zmax,pnum)
      real*8 gpressuremx(zmax,pnum,fnum)
      real*8 gcountmx(zmax,pnum,fnum)
      real*8 gprey(pnum)! sum of grazing rate*grazer biomass from each functional group
c      real*8 gpreytemp(pnum)! grazing rate*grazer biomass for just func gp grazer
      real*8 ggrazer(zmax,fnum)
      real*8 ggrazerbio(zmax,pnum)
c      real*8 testing(zmax,fnum)
      real*8 Ctemp
      real*8 clearance(zmax,pnum,fnum)
      real*8 capture(zmax,pnum,fnum)!capture efficiency, probability (so btw. 0 and 1)
      real*8 capexp(fnum) ! exponents to decide shape of capture efficiency
 


c *** Variable that record changes with time
      real*8 pchange(pnum,tnum)!prey biomass change with time
      real*8 zchange(zmax,tnum,fnum)! grazer biomass change with time
      real*8 Cdischange(tnum)!dissolved carbon change with time
      real*8 gtimechange(zmax,tnum,fnum) ! indicates which prey size class was grazed
      real*8 gchange(zmax,tnum,fnum) ! grazing pressure change with time

       real*8 zchangetemp1(zmax,fnum) ! temporary variable for recording
       real*8 zchangetemp2(zmax,fnum)
       real*8 pchangetemp1(pnum)
       real*8 pchangetemp2(pnum)
       real*8 Cdischangetemp




c ** Variables for saving things
      integer reclen


