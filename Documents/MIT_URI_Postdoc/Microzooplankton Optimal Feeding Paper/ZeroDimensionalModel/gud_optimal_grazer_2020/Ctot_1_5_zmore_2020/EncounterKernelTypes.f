c ************* Encounter Kernel for Specific Grazers **************
c *** This subroutine calculates the encounter kernel (m^3/s)
c *** for each predator encountering each type of prey based on 
c *** the relative size of predator and prey, just the size of the prey,
c *** empirically, or the size of the predator and prey

       subroutine EncounterKernelTypes (zmax, pnum,fnum,znumf,
     &     Dp, Dz, pr, zrmx, zvmx, pi,minind,v,alpha,
     &     functype,
     &     beta,encounterflag)

c *** Definig variables
       implicit none
       integer zmax
       integer pnum
       integer fnum
       integer znumf(fnum)

       double precision pi
       real*8 Dp(pnum)
       real*8 Dz(zmax,fnum)
       real*8 pr(pnum)
       real*8 zrmx(zmax,fnum)
       real*8 zvmx(zmax,fnum)
       real*8 minind(fnum)
       real*8 v(zmax,fnum)
       real*8 alpha
       integer functype(fnum)
       real*8 beta(zmax,pnum,fnum)
       integer encounterflag

       real*8 betatemp
       real*8 betatempv
       integer ind1, ind2
       integer l,j,k

c        print*,'Entering encounter kernel'

       do 700 l=1,fnum
             ind1 = minind(l)
c          print*,'ind1',ind1
       do 5000 j =1, znumf(l) 
              ind2 = ind1 +j - 1
c       print*,'zrmx(ind2,l)*10^6',zrmx(ind2,l)*10.0**6.0
          do 600 k = 1,pnum

c   *********** Encounter Kernel from swimming **************


c ********* Encounter kernel based on prey size (i.e., direct interception)
       if (encounterflag .eq. 1) then ! based just on prey size
         betatempv=(3.0/2.0)*pi*(pr(k)**2.0)
     &         *(v(ind2,l)**2.0+0.0)**(1.0/2.0)

c        betatempv=pi*((pr(k)+zrmx(ind2,l))**2.0)! based on prey+ predator
c     &         *(v(ind2,l)**2.0+0.0)**(1.0/2.0)
c


c ********  Piecewise equation with smooth transition based on relative size
c **** of predator and prey
      else if (encounterflag .eq. 2) then
c           Begining piecewise equation with smooth transition
             if (pr(k) .lt. zrmx(ind2,l)*alpha) then
         betatempv=(3.0/2.0)*pi*(pr(k)**2.0)
     &      *(v(ind2,l)**2.0+0.0)**(1.0/2.0)

         else if (pr(k).ge.zrmx(ind2,l)*alpha
     &     .and.pr(k).lt.zrmx(ind2,l)) then

           betatempv = (1-pr(k)/zrmx(ind2,l)) / (1-alpha) 
     &          *pi*(1.5)*v(ind2,l)*(pr(k)**2.0) +
     &       (1-((1-pr(k)/zrmx(ind2,l))/(1-alpha)))
     &       * pi*v(ind2,l)*((pr(k) + zrmx(ind2,l))**2.0)

             else if (pr(k) .ge. zrmx(ind2,l)) then
           betatempv = pi*v(ind2,l)*(pr(k) + zrmx(ind2,l))**2.0

c             else if (pr(k) .ge. zrmx(ind2,l)*alpha) then
c           betatempv = pi*v(ind2,l)*(pr(k) + zrmx(ind2,l))**2.0
c

             end if!if encounterflag = 2 and using varying beta


c ******* Using derivation cited in Kiorboe and Titelman (1998, pg. 1617)
 
      else if (encounterflag .eq. 3) then 
       betatempv = pi*v(ind2,l)
     &           *(pr(k)**2.0)*(3.0*zrmx(ind2,l)+2.0*pr(k))
     &            /(2.0*(pr(k)+zrmx(ind2,l)))




c ******* Empirical derivations
       else if (encounterflag .eq. 4) then !Empirical derivation

        if (functype(l) .eq. 1) then !dinoflagellate, (from Jeong et al., 2010)
           betatempv = (10.0**((4.18+log10(zvmx(ind2,l)*10.0**18.0)
     &  *(-1.08))))*(10.0**5.0)/3600.0*zvmx(ind2,l)

        elseif (functype(l) .eq. 2) then !ciliate (from Hansen et al, 1997)
          betatempv = (10.0**5.32)/(10.0**18.0)/3600.0
     &                 *(zvmx(ind2,l)*10.0**18.0)**0.98


        elseif (functype(l) .eq. 3) then !generic microzooplankter
        betatempv = (10.0**(-12.27))/3600.0
     &   *((4.0/3.0*pi*((zrmx(ind2,l)*10.0**6.0)**3.0))**0.82)

c        print*,'empirical relationshiop'
        end if! for empirical relationships



c ********* Based just on size of predator+prey all the time ************
        else if (encounterflag .eq. 5) then 
           betatempv = pi*v(ind2,l)*(pr(k) + zrmx(ind2,l))**2.0



      endif!encounterflag
c        print*,'betatempv',betatempv

c          *** Encounter kerenel from diffusion ***
           betatemp = 4.0*pi*(Dp(k)+Dz(ind2,l) )*(pr(k) + zrmx(ind2,l))

c          *** Adding together two encounter kernels ***
            beta(ind2,k,l) = betatempv + betatemp

c        print*,'beta from diffusion in loops',betatemp
c        print*,'beta from swimming in loops',betatempv

600       continue
5000      continue 
700     continue

c        print*,'betatemp',betatemp
c       print*, 'beta = ', beta
c       print*, 'Leaving encounter kernel'
c        print*,'beta from diffusion',betatemp
c        print*,'beta from swimming',betatempv
       
      end ! EncounterKernelTypes subroutine




