c ************  Grazing based on clearance rate **************

c *** Use this code to determine which size class of prey 
c *** is grazed by each predator based on the clearance 
c *** rate scaled by the phytoplankotn (prey) biomass


      subroutine GrazingRateChoiceIV(pnum,zmax,fnum,
     &  minind,maxind,
     &  pCcell, zCcell, 
     &  pCsize,!carbon concentration for each size class
     &  zCsize, !the carbon concentration for each size class
     &  beta, h, gmax, clearance,
     &  capture, grazeflag,prefpreyclass,
     &  gtimetemp,gtemp,!gcount
     &   gprey,ggrazer)!testing)

c *** This subroutine calculates the preferred grazing rate 
c *** based on the clearance rate.  It's based off the code
c *** "GrazingRateIII.f"

c *** Defining parameters
      implicit none


c *** i, j, initz, and initp are all used for loops
c    maxi and maxj contain indices for the matrices that keep 
c    track of the grazing rate and the prey size class grazed

      integer pnum
      integer zmax
      integer fnum
      real*8 minind(fnum)
      real*8 maxind(fnum)

      real*8 pCcell(pnum)
      real*8 zCcell(zmax,fnum)
      real*8 pCsize(pnum)
      real*8 zCsize(zmax,fnum)
      real*8 beta(zmax,pnum,fnum)
      real*8 h(zmax,pnum,fnum)
      real*8 gmax(zmax,pnum,fnum)
      real*8 clearance(zmax,pnum,fnum)
      real*8 capture(zmax,pnum,fnum)!capture probability
      real*8 gpressure(zmax,pnum)! grazing rate
      integer prefpreyclass(fnum)! preferred prey class--only imp. if grazeflag=3
c      real*8 gcount(zmax,pnum)
      real*8 gtimetemp(zmax,fnum)
      real*8 ggrazerbio(zmax,pnum)! grazing rate*grazer biomass
      real*8 ggrazer(zmax,fnum) ! grazing rate*grazer biomass FOR EACH GRAZER
      real*8 gprey(pnum) ! grazing rate*grazer biomass for each prey
      real*8 gtemp(zmax,pnum)
      real*8 tempval!temporary varible in code for finding max grazing
      integer maxi(zmax)
      integer maxj(zmax)


c *** Flags
      integer grazeflag! to determine how preferred grazed size class is chosen 
      integer grazeflag1! have to have these extra flags because the flag
c            value changed during the subroutine--weird
      integer grazeflag2
 

c *** Indices for loops 
      integer ind1, ind2
      integer indi, indj
      integer i,j,jj,ii,ij,f

      

c       print*, 'Entering grazing rate subroutine'

c       print*,'grazeflag in subroutine',grazeflag

c *** Defining flags
       grazeflag1 = grazeflag
       grazeflag2 = grazeflag

         

******************************************************
c *** Loop to go through each functional group ******
******************************************************

         do 220 f = 1,fnum
          ind1 = minind(f) ! Defining indices for each functional group
          ind2 = maxind(f)
        
c          print*,'f',f
c          print*,'ind1',ind1
c          print*,'ind2',ind2
c
c          print*,' '
c          print*,' '
c

c *** Initializing things to zero

      do 900 i = 1,zmax
c            gtimetemp(i,f) = 0.0
        do 1000 j = 1,pnum
            gpressure(i,j) = 0.0
            ggrazerbio(i,j) = 0.0
            gtemp(i,j) = 0.0
1000    continue
900   continue

c            print*,'gtemp',gtemp

c   *** If statement determining how size class grazed is determined
c       print*,'clearance(:,:,1)',clearance(:,:,1)
c       print*,'clearance(:,:,2)',clearance(:,:,2)

c        print*,'pCsize',pCsize
 
   

      if (grazeflag1 .eq. 1) then ! size class based on actual highest intake
      do 700 i=ind1,ind2
       do 800 j=1,pnum
         gtemp(i,j) = gmax(i,j,f)*(pCsize(j)/(pCsize(j) !<--Changed in 2020, too
c     &      + zCcell(i,f)/beta(i,j,f)/h(i,j,f)/capture(i,j,f) ))
     &      + pCcell(j)/beta(i,j,f)/h(i,j,f)/capture(i,j,f) ))!<- This is what changed in 2020
800    continue
700   continue

       elseif (grazeflag1 .eq. 2) then ! size class based on scaled clearance rate
c          print*,'grazed size class chosen by clearance'
      do 750 i = ind1,ind2
        do 850 j = 1,pnum
         gtemp(i,j) = clearance(i,j,f)*pCsize(j)*capture(i,j,f)
850   continue
750   continue


       elseif (grazeflag1 .eq. 3) then!pre-selected preferred size class
        do 775 i = ind1,ind2
           do 875 j = 1,pnum
              if (j .eq. prefpreyclass(f)) then
         gtemp(i,j) = gmax(i,j,f)*(pCsize(j)/(pCsize(j) !<--Changed in 2020, too
c     &      + zCcell(i,f)/beta(i,j,f)/h(i,j,f)/capture(i,j,f) ))
     &      + pCcell(j)/beta(i,j,f)/h(i,j,f)/capture(i,j,f) ))!<- This is what changed in 2020
              else
         gtemp(i,j) = 0.0
             end if
875        continue
775     continue


       end if
    
c          print*,'grazeflag1 after first if',grazeflag1
c          print*,'grazeflag after first if',grazeflag
c
c         print*, 'gtemp(1,:)',gtemp(1,:)
c         print*,'gtemp(2,:)', gtemp(2,:)
c         print*,'gtemp(3,:)', gtemp(3,:)
c         print*,'gtemp(4,:)', gtemp(4,:)
c         print*,'gtemp(5,:)', gtemp(5,:)
c         print*,'gtemp(6,:)', gtemp(6,:)
c         print*,'gtemp(7,:)', gtemp(7,:)
c         print*,'gtemp(8,:)', gtemp(8,:)
c         print*,'gtemp(9,:)', gtemp(9,:)
c


c       print*,'initial gpressure',gpressure
c        print*,'ggrazerbio initial',ggrazerbio
   
c          print*,'grazeflag1 after intializing to zero',grazeflag1


c *** Finding the value and indices for the max grazing rate
         tempval = gtemp(1,1) !setting initial max value
         maxi = 0
         maxj = 0

c         print*, 'maxi initial',maxi
c        print*, 'maxj initial',maxj

      do 1100 i = ind1,ind2
        do 1200 j = 1,pnum
          if (j .eq. 1) then
               tempval = gtemp(i,1)
               maxj(i) = j
         end if
            maxi(i) = i
          if ( gtemp(i,j) .gt. tempval )then
             maxj(i) = j
             tempval = gtemp(i,j)
          end if
c                print*,'tempval',tempval

1200    continue
1100    continue

c        print*,'maxi final',maxi
c        print*,'maxj final',maxj

c      print*,'gpressure(1,:)',gpressure(1,:)
c      print*,'gpressure(2,:)',gpressure(2,:)
c      print*,'gpressure(3,:)',gpressure(3,:)
c      print*,'gpressure(4,:)',gpressure(4,:)
c      print*,'gpressure(5,:)',gpressure(5,:)
c      print*,'gpressure(6,:)',gpressure(6,:)
c      print*,'gpressure(7,:)',gpressure(7,:)
c  
       
c *** Putting in the max grazing pressure and keeping track
c      of the size class for the appropriate grazer and prey 
c      size classes that lead to the max grazing rate


      if (grazeflag2 .eq. 1 .or. grazeflag2 .eq. 3) then! size class grazed det. by max biomass incr. from grazing or have a pre-selected preferred prey size class
c          print*,'grazeflag2 in if',grazeflag2


      do 1300 i= ind1,ind2
               indi = maxi(i)
               indj = maxj(i)

c         print*,'indi',indi
c         print*,'indj',indj
c         print*,'gtemp(indi,:)',gtemp(indi,:)

      gpressure(indi, indj) = gtemp(indi, indj)
c      gcount(indi, indj) = gcount(indi,indj) + 1.0
      gtimetemp(i,f) = indj ! I think this is the prey size class that each grazer consumes

c         print*,'gpressure(indi,:)',gpressure(indi,:)
c         print*,'gpressure',gpressure

1300    continue
        else if (grazeflag2 .eq. 2) then! size class grazed det. by clearance rate
c           print*,'grazeflag2',grazeflag2
      do 1350 i= ind1,ind2
               indi = maxi(i)
               indj = maxj(i)
       gpressure(indi, indj) = 
     &                1/h(indi,indj,f)
     &    *(pCsize(indj)/(pCsize(indj) 
c     &      + zCcell(indi,f)/beta(indi,indj,f)/h(indi,indj,f)
     &      + pCcell(indj)/beta(indi,indj,f)/h(indi,indj,f)!<- Changed in 2020
     &      /capture(indi,indj,f)))


c      gcount(indi, indj) = gcount(indi,indj) + 1.0
      gtimetemp(i,f) = indj ! I think this is the prey size class that each grazer consumes
1350    continue

      end if !grazeflag2

c          print*,'grazeflag2 after second if',grazeflag2
c          print*,'grazeflag after second if',grazeflag

 
c      print*,'gtimetemp',gtimetemp


c      print*, 'gpressure =', gpressure


c      print*,'gpressure(1,:)',gpressure(1,:)
c      print*,'gpressure(2,:)',gpressure(2,:)
c      print*,'gpressure(3,:)',gpressure(3,:)
c      print*,'gpressure(4,:)',gpressure(4,:)
c      print*,'gpressure(5,:)',gpressure(5,:)
c      print*,'gpressure(6,:)',gpressure(6,:)
c      print*,'gpressure(7,:)',gpressure(7,:)
c      print*,'gpressure(8,:)',gpressure(8,:)
c      print*,'gpressure(9,:)',gpressure(9,:)
c
  
c       print*,'gcount',gcount
c       print*, 'gtemp =', gtemp
cc       print*, 'temmpval = ', tempval
c       print*, 'maxi =', maxi, 'maxj=', maxj
c       print*,' gcount(indi,indj)', gcount(indi,indj)

c       print*,'zCcell',zCcell



c *** Mutliplying grazing rate by grazer biomass because
c * later the values are summed together to get grazing pressure
c * on prey and by grazers, which changes the dimensions.  So,
c * it's just easier to do the multiplying by grazer biomass now

c         print*,'gpressure before loop',gpressure
          do 350 i= ind1,ind2
              do 450 j=1,pnum
c                  print*,'i',i
c                  print*,'gpressure(i,:) initial',gpressure(i,:)
c                  print*,'zCsize(i)',zCsize(i)


                  ggrazerbio(i,j) = gpressure(i,j)*zCsize(i,f)


c                  print*,'ggrazerbio(i,:)',ggrazerbio(i,:)
450           continue
350       continue


c      print*,'gpressure(1,:)',gpressure(1,:)
c      print*,'gpressure(2,:)',gpressure(2,:)
c      print*,'gpressure(3,:)',gpressure(3,:)
c      print*,'gpressure(4,:)',gpressure(4,:)
c       print*,'gpressure(5,:)',gpressure(5,:)
c       print*,'gpressure(6,:)',gpressure(6,:)
c      print*,'gpressure(7,:)',gpressure(7,:)
c         
          
c          print*,'grazer biomass',zCcell

c         print*,' '
c         print*,'ggrazerbio(1,:)',ggrazerbio(1,:)
c         print*,'ggrazerbio(2,:)',ggrazerbio(2,:)
c         print*,'ggrazerbio(3,:)',ggrazerbio(3,:)
c         print*,'ggrazerbio(4,:)',ggrazerbio(4,:)
c          print*,'ggrazerbio(5,:)',ggrazerbio(5,:)
c          print*,'ggrazerbio(6,:)',ggrazerbio(6,:)
c         print*,'ggrazerbio(7,:)',ggrazerbio(7,:)
c          print*,'ggrazerbio(8,:)',ggrazerbio(8,:)
c         print*,'ggrazerbio(9,:)',ggrazerbio(9,:)
c
c      

c * Finding the grazing pressure on each prey size class 
c * and the grazing pressure induced by each grazer size class

c           print*,'ggrazerbio',ggrazerbio

c           print*,'pnum',pnum

          do 150 j=1,pnum
           do 155 jj = ind1,ind2
c          print*,'j',j
c          print*,'ggrazerbio(:,j)',ggrazerbio(:,j)

c          gprey(j) = gprey(j) + sum(ggrazerbio(ind1:ind2,j))

          gprey(j) = gprey(j) + ggrazerbio(jj,j)


155          continue

c          print*,'gprey(j)',gprey(j)

150       continue
   


c          print*,'znum',znum

          do 250 ii= ind1,ind2
            do 255 ij = 1,pnum
c             print*,'ii',ii
c             print*,'ggrazerbio(ii,:) before row sum',ggrazerbio(ii,:)
c              ggrazer(ii,f) = sum(ggrazerbio(ii,1:pnum))
              ggrazer(ii,f) = ggrazer(ii,f) + ggrazerbio(ii,ij)

c             print*,'ggrazerbio(ii,:) after  row sum',ggrazerbio(ii,:)

c               print*,'ggrazer(ii,f)',ggrazer(ii,f)
255          continue
250       continue

c          print*,'grazeflag1 at end of subroutine',grazeflag1
c          print*,'grazeflag2 at end of subroutine',grazeflag2


          grazeflag = grazeflag1


c          print*,'grazeflag at end of subroutine',grazeflag
c         print*,'ggrazer',ggrazer
c

c         print*,'grazing on prey',gprey
c         print*,'grazing by grazer',ggrazer


c
220     continue! do loop going through all functional groups

c       print*, 'Leaving grazing rate subroutine'

       end ! grazing_rate subroutine     







