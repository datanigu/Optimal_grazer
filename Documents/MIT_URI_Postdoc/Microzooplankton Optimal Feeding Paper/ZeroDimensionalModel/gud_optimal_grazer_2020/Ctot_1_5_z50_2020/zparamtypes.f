c ******** Grazer physiological and associated parameters ************

        subroutine zparamtypes(fnum,zmax,znumf,zrmx,zvmx,pi,pnum,pr,
     &   minind, maxind,tmp,etatmp,
     &   gtmpcoeff, gmaxexp,gE,kBoltz,
     &   bodylengths,functype,
     &   efficiency,zm0,RE,zCcell,
     &   v,Dz,gmax,zm,R0, Rv,R,swimregflag)
c    &   betav)

c Declaring variables
       implicit none
       integer fnum ! number of functional groups
       integer zmax !absolute max number of size classes
       integer znumf(fnum) !number of size classes in each functional group
       real*8 zrmx(zmax,fnum)! radius, m
       real*8 zvmx(zmax,fnum)! volume, m^3
       double precision pi !pi
       integer pnum !number of phytoplankton groups
       real*8 pr(pnum) ! radius of phytoplankton, m
       real*8 minind(fnum)!index of min size class out of all size classes
       real*8 maxind(fnum) ! max number of size classes out of all classes
       real*8 v(zmax,fnum)!swimming speeds
       real*8 Dz(zmax,fnum) !diffusion coefficient for grazers
       real*8 gmax(zmax,pnum,fnum) ! max grazing rate , 1/s
       real*8 gE !grazing activation energy, eV
       real*8 gmax0 !max grazing coefficient includng temperature-dep
       real*8 gtmpcoeff ! max grazing coeff not including temp-dep
       real*8 gmaxexp !max grazing exponent 
       real*8 kBoltz ! Boltzmann's contant = 8.617*10^-5 eV/K (where K = Kelvin)
       real*8 zm0(fnum) !basal mortality rate (i.e., not due to swimming), 1/s      
       real*8 zm(zmax,fnum)! 1/s, total grazing mortality, from swimming and basal 
       real*8 z0R(zmax,fnum) ! O2/cell/hr, basal respiration for grazer
       real*8 R0(zmax,fnum) !1/s, specific basal respiration for grazer  
       real*8 efficiency !efficiency of ciliary/flagellar propulsion
       real*8 Rv(zmax,fnum) ! 1/s, specific respiration from swimming
       real*8 RE ! eV, actiation energy for respiration
       real*8 zCcell(zmax,fnum)!carbon per cell, gC/cell
       real*8 tmp !temperature, C
       real*8 etatmp !temperature-dep kinematic viscosity for respiration, g/m/s, 
                    ! from Sharqawy et al., 2010, equation 22
       real*8 R(zmax,fnum) !1/s, total respiration (basal and swimming)  
       real*8 bodylengths(fnum) ! mutliple of bodylengths swim per second
c            only important if swimregflag not 1
       integer swimregflag(fnum) ! determine if use regression to find swim speed
       integer functype(fnum) ! determine what functional type (1=dino,2=ciliate,
c                      3 = generic microzoopl)

       integer q,st !indices for do loops
       integer ind1, ind2 !for use in indexing WITHIN loop
       integer di, ci, gi, gp_i, dp_i, cp_i !for loops



c ***** Begin loop to go through functional groups
      do 500 q = 1,fnum 
         ind1 = minind(q)
         ind2 = maxind(q) 
c         print*,'starting index',ind1
c         print*,'ending index',ind2


c ************** Dinoflagellate ***********************
        if (functype(q) .eq. 1) then !dino
        do 501 di= ind1,ind2
         do 5011 dp_i = 1,pnum

c *** grazing
          gmax(di,dp_i,q) = (0.178*(2.0*zrmx(di,q)*10.0**6.0)+2.89)
     &          *(10.0**(-9.0))/24.0/3600.0/zCcell(di,q)
c *** Swimming speed
              if (swimregflag(q) .eq. 1) then !use regression equation
          v(di,q) = 1.1879*((zvmx(di,q)*10.0**18.0)**0.5094)
     &        *(10.0**(-6.0))
              else !swim speed based on multiple of body lengths
          v(di,q) = zrmx(di,q)*2.0*bodylengths(q)
              end if !swimregflag
5011    continue
501     continue


c ************** Ciliate *******************
       elseif (functype(q) .eq. 2) then !ciliate
        do 502 ci = ind1,ind2
         do 5021 cp_i = 1,pnum
c *** grazing
         gmax(ci,cp_i,q)=(10.0**0.1)/3600.0*
     &      (zvmx(ci,q)*10.0**18)**(-0.2)
c **** swimming
              if (swimregflag(q) .eq. 1) then
         v(ci,q) = (127.1603*(zvmx(ci,q)*10.0**18.0)**0.1656)
     &            *(10.0**(-6.0))
              else
          v(ci,q) = zrmx(ci,q)*2.0*bodylengths(q)
              end if !swimregflag
5021     continue
502      continue



c ***************** Generic microzooplankton ************
        elseif (functype(q). eq. 3) then !generic microzooplankton
           do 503 gi = ind1,ind2
            do 504 gp_i = 1,pnum
c *** Grazing
c          gmax(gi,q) = gtmpcoeff*((2.0*zrmx(gi,q)*(10**6))**gmaxexp)
c     &   *exp(-gE/kBoltz/(tmp+273.15)) 
c     &   /24.0/3600.0
         gmax(gi,gp_i,q) = 0.0035*((zrmx(gi,q)/pr(gp_i))** (-1.9431))
504       continue
c *** Swimming
          v(gi,q) = zrmx(gi,q)*2.0*bodylengths(q)
503       continue 
        endif !functype


c ***** Parameters that have same equation regardless of func type
        do 700 st = ind1,ind2 ! going through each size class
c *** Diffusion coefficient
        Dz(st,q) = ((tmp+273.15)*(1.38e-20))/(6*pi*etatmp*zrmx(st,q))

c *** Mortality and respiration 
      zm(st,q) = zm0(q) + zm0(q)*v(st,q)
c
      z0R(st,q)=10.0**(0.75*log(zvmx(st,q)*10.0**18.0)/log(10.0)-4.09)

c       print*,'z0R',z0R

c      R0(st,q) = z0R(st,q)/(zvmx(st,q)*10.0**18)*(1./10.0**9)*
c     &      (1./0.08206/(tmp+273.15))*
c     &      12.0*(10.**12)/0.216/3600.0!1/sec, basal respiration rate--this is WRONG



c      R0(st,q) = z0R(st,q)/(10.0**9.0)*(1.0/0.08206/(tmp+273.15))*12.0*
c     &       (10.0**12.0)/(0.216*(zvmx(st,q)*10.0**18.0)**0.939)/3600.0! CORRECT! Based on Fenchel and Finlay (1983)

      R0 = 0.0 ! Use if respiration scales and is accounted for in ingestion (from Verity 1985)

c      Rv(st,q) = ((3.0*pi*zrmx(st,q)*2.0*(v(st,q)**2.0))*etatmp/1000.0
c     &        *(10.0**(-3.0))/20.2/(tmp+273.15) 
c     &        / 0.082*12.0/zCcell(st,q)/efficiency)
ccc     &         *exp(-RE/kBoltz/(273.15+tmp))
c

      Rv(st,q) = ((3.0*pi*zrmx(st,q)*2.0*(v(st,q)**2.0))*etatmp/1000.0
     &        *(10.0**(-3.0))/20.2/(tmp+273.15) 
     &        / 0.082*12.0/zCcell(st,q)/efficiency)

c
       R(st,q) = R0(st,q) + Rv(st,q)

      
700    continue ! going through each size class


500    continue ! goingthrough each functional group


c        print*,'variables in subroutine'
c       print*, 'v=',v
c       print*, 'betav',betav
c       print*, 'Dz=',Dz
c        print*,'gmax col 1=',gmax(:,1)
c        print*,'gmax col 2=',gmax(:,2)
c        print*, 'zm=',zm
c        print*,'z0R',z0R
c        print*,'R0 in subroutine',R0
c        print*,'zrmx in subroutine',zrmx
c        print*,'v',v
c        print*,'etatmp',etatmp
c        print*,'Rv in subroutine',Rv
c        print*, 'R=',R
c
c        print*,'leaving zparamtypes.f subroutine'

c       return
       end
