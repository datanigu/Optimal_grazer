c *********** Handling time
c ** This subroutine calculates the handling time, h, for 
c ** each grazer consuming each different typ of prey

      subroutine handlingIII(zmax, pnum,fnum,znumf,
     &   gmax, pCcell, zCcell,zrmx,pr,
     &    minind,maxind, hflag,
     &     h)

c *** Defining variables
      implicit none
      integer zmax, pnum,fnum, hflag
      integer znumf(fnum)
      real*8 minind(fnum)
      real*8 maxind(fnum)
       
      real*8 gmax(zmax,fnum)
      real*8 zCcell(zmax,fnum)
      real*8 pCcell(pnum)
      real*8 zrmx(zmax,fnum)
      real*8 pr(pnum)
      real*8 h(zmax,pnum,fnum)

      integer  k, j,l,kk,jj,ll,lll,jjj,kkkk!for do loops
      integer ind1, ind2 !integers for indexing within do loops


c         print*, 'Entering handling time'
c         print*, 'gmax',gmax
c         print*, 'pCcell',pCcell
c         print*, 'zCcell',zCcell
c          print*,'hflag', hflag
c
c *** Calculating handling time
       if (hflag .eq. 1) then ! handling based on max grazing rate
      do 200 l = 1,fnum !each functional group
        ind1 = minind(l)
c         print*,'minind(l)',minind(l)
c         print*,'ind1',ind1
      do 300 j = 1, znumf(l) !each row is different predator
          ind2 = ind1 + j -1
c            print*,'ind2',ind2
         do 400 k = 1, pnum ! each column is different prey
          h(ind2,k,l) = 1/gmax(ind2,l)*( pCcell(k)/zCcell(ind2,l) )  
400      continue 
300   continue 
200   continue

           elseif (hflag .eq. 2) then ! handling based on other data and has threshold
      do 250 ll = 1,fnum !each functional group
        ind1 = minind(ll)
      do 350 jj = 1, znumf(ll) !each row is different predator
          ind2 = ind1 + jj -1
         do 450 kk = 1, pnum ! each column is different prey
         h(ind2,kk,ll) = 2.16e5*(zrmx(ind2,ll)/pr(kk))**(-1.7782) + 0.5
450      continue 
350   continue 
250   continue


           elseif (hflag .eq. 3) then ! handling based on assumptions and NO threshold
      do 275 lll = 1,fnum !each functional group
        ind1 = minind(lll)
      do 375 jjj = 1, znumf(lll) !each row is different predator
          ind2 = ind1 + jjj -1
         do 475 kkkk = 1, pnum ! each column is different prey
         h(ind2,kkkk,lll) = 382.981*(zrmx(ind2,lll)/pr(kkkk))**(-1.006)
475      continue 
375   continue 
275   continue


         endif
 
c
c       print*,'h(1,1,1)', h(1,1,1)
c       print*,'zrmx(1,1)/pr(1)',zrmx(1,1)/pr(1)
c       print*,'h2(1,1,1) in subroutine',
c     &      2.16e5*(zrmx(1,1)/pr(1))**(-1.7782) + 0.5
c      print*, 'handling time group 1 =', h(:,:,1)
c      print*, 'handling time group 2 =', h(:,:,2)
c      print*,' '
c      print*,' '


c      print*, 'Leaving handling time'
c         print*, 'gmax',gmax
c         print*, 'pCcell',pCcell
c           print*, 'zCcell',zCcell


      end! handling subroutine



