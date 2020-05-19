% Debugging zero-d optimal grazer model
mainFolder = '/Users/dtaniguchi/Documents/MIT_URI_Postdoc/Microzooplankton Optimal Feeding Paper/ZeroDimensionalModel/gud_optimal_grazer_2020';
subdirs_order = {'Ctot_1_5_z05_2020','Ctot_1_5_z50_2020','Ctot_1_5_z200_2020','Ctot_1_5_zmany_200_2020'};

cd(mainFolder)
    
    kk = 2;
    %% Name of this subdirectory
subdname = subdirs_order{kk};
eval(['addpath ', subdname])

%% Reading in binary files
    % pnum = number of prey size classes
    fid=fopen('pnum.bin','r','n');
    pnum = fread(fid,'integer*8');
    fclose(fid);
    
    % zmax = max number of grazer size classes
    fid=fopen('zmax.bin','r','n');
    zmax = fread(fid,'integer*8');
    fclose(fid);
    
    % fnum = number of grazer functional groups
    fid=fopen('fnum.bin','r','n');
    fnum = fread(fid,'integer*8');
    fclose(fid);
    
    % dt = time step between each time point, s
    fid=fopen('dt.bin','r','n');
    dt = fread(fid,'real*8');
    fclose(fid);
    
    % tnum = the time step values, s
    fid=fopen('tnum.bin','r','n');
    tnum= fread(fid,'integer*8');
    fclose(fid);  
    tnum_all = tnum; 
    
savestep = 1;% will get written over if not every time step saved

    %if all_steps_flag ~=1        
        % tcount = number of time steps actually recorded
        fid=fopen('tcount.bin','r','n');
        tcount= fread(fid,'integer*8');
        fclose(fid);
        tnum = tcount;
        
        % savestep = number of time steps between each saved value
        fid=fopen('savestep.bin','r','n');
        savestep= fread(fid,'integer*8');
        fclose(fid);     

    %end %all_steps_flag
    
    % pr = prey radius values
    fid=fopen('pr.bin','r','n');
    pr= fread(fid,'real*8');
    fclose(fid);
    
    % zrmx = predator radius values
    fid = fopen('zrmx.bin','r','n');
    zrmxtemp=fread(fid,'real*8');
    zrmx = reshape(zrmxtemp,zmax,fnum);
    fclose(fid);
    
    % zsizewidth = log based scale used to determine predator radius values
        fid=fopen('zsizewidth.bin','r','n');
        zsizewidth= fread(fid,'real*8');
        fclose(fid);
     
    % minind = vector with indices for min size class for each predator
    % functional group
    fid = fopen('minind.bin','r','n');
    minind = fread(fid,'real*8');
    fclose(fid);
    
    %maxind = vector with indices for max size class for each predator
    % in each functional group
    fid = fopen('maxind.bin','r','n');
    maxind = fread(fid,'real*8');
    fclose(fid);
    
    %pchange = change in prey biomass with time, gC/m^3
    % where each row is a time step, and each column is a prey size class
    fid=fopen('pchange.bin','r','n');
    pchangetemp = fread(fid,'real*8');
    fclose(fid);
    pchange = reshape(pchangetemp,pnum,tnum);
    
    % zchange = change in grazer biomass with time, gC/m^3
    % where each row is a time step, and each column is a prey size class,
    % and the third dimension is for each size class
    fid=fopen('zchange.bin','r','n');
    zchangetemp = fread(fid,'real*8');
    fclose(fid);
    zchange = reshape(zchangetemp,zmax,tnum,fnum);
    
    
    % Cdischange = change in dissolved carbon biomass with time, gC/m^3
    % where each row is a time step, and each column is a prey size class
    fid=fopen('Cdischange.bin','r','n');
    Cdischange = fread(fid,'real*8');
    fclose(fid);
    
    % Ctot = total carbon = z_change(1,:) + p_change(1,:) + Cdis_change(1,:)
    % in units gCm^3
    fid=fopen('Ctot.bin','r','n');
    Ctot = fread(fid,'real*8');
    fclose(fid);
    
    % v = swimming speed (m/s)
    fid=fopen('v.bin','r','n');
    v = fread(fid,'real*8');
    fclose(fid);

    % gtimechange = list of times each size class chosen by each grazer to
    % consume at each time step
    fid = fopen('gtimechange.bin','r','n');
    gtimechangetemp = fread(fid, 'real*8');
    fclose (fid);
    %if all_steps_flag ~=1
        gtimechange=reshape(gtimechangetemp,zmax,tnum_all,fnum);
%     else
%         gtimechange=reshape(gtimechangetemp,zmax,tnum,fnum);
%     end% all_steps_flag
    
    % grazeflag = if 1, then grazing based on Michaelis-Menten Derivation
    % of max grazing rate; otherwise based on clearance rate
    fid = fopen('grazeflag.bin','r','n');
    grazeflag = fread(fid, 'integer*8');
    fclose (fid);
    
    % encounterflag = 1 means encounter rate based solely on prey;
    % otherwise, encounter based on relative size of prey and grazer
    fid = fopen('encounterflag.bin','r','n');
    encounterflag = fread(fid, 'integer*8');
    fclose (fid);
 
    
    % captureflag = 1 means 100% capture efficiency probability for grazers consuming all sized prey; 
    % otherwise, capture probability decreases
    % with increasing prey size
    fid = fopen('captureflag.bin','r','n');
    captureflag = fread(fid, 'integer*8');
    fclose (fid);
    
    % mortalityflag = 1 means linear grazer mortality, 2 means quadratic grazer mortality; 
    fid = fopen('mortalityflag.bin','r','n');
    mortalityflag = fread(fid, 'integer*8');
    fclose (fid);

     % muflag = 1 means linear change in phyto growth with size; 2 means unimodal growth
    fid = fopen('muflag.bin','r','n');
    muflag = fread(fid, 'integer*8');
    fclose (fid);

    %swimflag = 1 means constant swim speed; 2 means changes with grazing intake
    fid = fopen('swimflag.bin','r','n');
    swimflag = fread(fid, 'integer*8');
    fclose (fid);
    
    %swimregflag = 1 means based on regression; 2 means based on multiple of cell DIAMETER (NOT radius)
    fid = fopen('swimregflag.bin','r','n');
    swimregflag = fread(fid, 'integer*8');
    fclose (fid);
    
    % Carbon content per phytoplankton cell, g/cell
    fid = fopen('pCcell.bin','r','n');
    pCcell = fread(fid,'real*8');
    fclose(fid);
    
    
    % Carbon concent per grazer cell, g/cell
    fid = fopen('zCcell.bin','r','n');
    zCcell = fread(fid,'real*8');
    fclose(fid);  
    
    
 %% Printing values of interest
 
 disp('pCcell is ')
 disp(pCcell)
 disp('zCcell is ')
 disp(zCcell)
    
    
 
%% Calculating values
%pradius = 0.50001096661654242*10^-6;
pradius = [0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256];%um 
pvol = 4/3*pi.*(pradius.^3);%um^3
%zradius = 2.50*10^-6;
zradius = [2.5, 25, 100];%um
zvol = 4/3*pi.*(zradius.^3);%um^3

Qz = 0.216*(zvol).^0.939;%pg C/cell
Qz = Qz.*(10^-12)/12*1000;%mmolC/cell
Qz_mx = repmat(Qz,length(pradius),1);

Qp = 0.216*(pvol).^0.939;%pg C/cell
Qp = Qp*10^-12/12*1000;%mmolC/cell
Qp_mx = repmat(Qp,length(zradius),1)';

[x,y] = meshgrid(zradius, pradius);
handling = 3600.*(x./y).^(-1.27);%s
%phandling = 3600.*(zradius./pradius).^-1.27;

captureprob = (1-y./(y+x)).^3;%dimensionless
% captureprob = (1-pradius/(pradius+zradius))^3;

%bodylength = 20;%number of body lengths swim/second
bodylength = 10;%number of body lengths swim/second
eta = 1.1385e-6;%dynamic viscosity, g/um/s, taken from gud_tempfunc.F, which
    % should hopefully match etatmp below, but 10^6 times larger (note units)

vswim = 2*bodylength*zradius;%um/s
vswim_mx = repmat(vswim,length(pradius),1);
betav=3.1419527.*vswim_mx.*(x+y).^2;%um^3/s
%betav=3.1419527*vswim*(pradius+zradius)^2;%um^3/s
%betav = betav*10^-18;%m^3/s

sal = 0.00035;
tmp = 15;
etaW =  4.2844e-5+(0.157*(tmp+64.993)^2.0-91.296)^(-1.0);
etatmp = etaW* ...
      ( 1+ (1.541+(1.998*10^(-2.0))*tmp-(9.52*10^(-5.0))*tmp^2.0)* ...
          (sal) + ...
          ((7.974-7.561*10^(-2.0))*tmp+...
              (4.724*10^(-4.0))*tmp^2.0)*(sal^2.0))*1000;%g/m/s
          
diffcoeff_z = ((tmp+273.15)*(1.38e-20))./(6*pi*etatmp*zradius*10^-6);
diffcoeff_zmx = repmat(diffcoeff_z,length(pradius),1);
diffcoeff_p = ((tmp+273.15)*(1.38e-20))./(6.0*pi*etatmp*pradius*10^-6);
diffcoeff_pmx = repmat(diffcoeff_p,length(zradius),1)';
betad_m3s = 4.0*3.1415927*(diffcoeff_zmx+diffcoeff_pmx ).*(y.*10^-6 + x.*10^-6);%m^3/s
betad = betad_m3s.*((10^-6)^3);%um^3/s

kz = Qz_mx./captureprob./handling./(betav+betad);%mmolC/um^3 
kz = kz.*(10^18);%umol C/L
kz = kz.*16./106;%umol N/L

kz_p = Qp_mx./captureprob./handling./(betav+betad);%mmolC/um^3
kz_p = kz_p.*(10^18);%umol C/L
kz_p = kz_p.*16./106;%umol N/L

% Other parameters, for completeness
mur = 1.36.*(2.*pradius).^-0.16;
muv = 1.3139.*pvol.^-0.0533;
ksr = 0.33.*(2.*pradius).^0.48;%nutrient half sat based on radius
ksv = 0.366.*pvol.^0.16;%nutrient half sat based on vol

%% Comparing what kz (grazing half saturation constant) looks like when use 
% carbon content of grazers compared to carbon content of prey

disp('kz');kz
disp('kz with Qp');kz_p

disp('max absolute kz-kz_p')
max(abs(kz-kz_p))

disp('max absolute (kz-kz_p)/kz_p')
max_diff = max(abs(kz-kz_p)./kz_p)

%% Based on min and max values of direct measurements of handling time from Boegnik and Arndt (200)
% handling time ranges from 3.87 s to 94.5 s

hmin = 3.87; %seconds
hmax = 94.5; % seconds

kz_hmin = Qp_mx./captureprob./hmin./(betav+betad)*1000*(10^18)/1000*16/106;%umol N/L
kz_hmax = Qp_mx./captureprob./hmax./(betav+betad)*1000*(10^18)/1000*16/106;%umol N/L

%% Grazing half saturation constant based on direct measurements of handling time
% in Handling_Time_Parameterization_II.m
% For direction measurements of handling time, coeff_direct_hand  =  382.9811, exp_direct_hand =  -0.3354
% That is, handling time in seconds = 382.9811 * (predvol/preyvol)^-0.3354

[x,y] = meshgrid(zradius, pradius);
zvol = 4/3*pi.*x.^3;
pvol = 4/3*pi.*y.^3;
h_direct = 382.9811.*(zvol./pvol).^(-0.3354);%s

kz_h_direct = Qp_mx./captureprob./h_direct./(betav+betad)*1000*(10^18)/1000*16/106;%umol N/L


%% Grazing half saturation constant based on inverse of max grazing rate, from relationship
% in Handling_Time_Parameterization_II.m
% For inverse of max grazing, handling time in seconds = 410.3568 * (predvol/preyvol)^0.5551
h_max_graze = 410.3568.*(zvol./pvol).^0.5551;
kz_max_graze = Qp_mx./captureprob./h_max_graze./(betav+betad)*1000*(10^18)/1000*16/106;%umol N/L



%% Grazing half saturation constant based on combo of direct measures of handling time and inverse of max grazing rate
% in Handling_Time_Parameterization_II.m
% For all data combined, handling time in seconds = 293.0759*(predvol/preyvol)^0.4368 
h_all = 293.0759.*(zvol./pvol).^0.4368;
kz_h_all = Qp_mx./captureprob./h_all./(betav+betad)*1000*(10^18)/1000*16/106;%umol N/L

%% Plotting grazing half saturation
figure(1);
current_axes=gca;
surf(kz_h_direct)
xlabel ('predator vol, \mum^3');ylabel ('prey vol, \mum^3');
zlabel ('kz, umol N/L');
title('kz based on direct measures of handling time');
set(current_axes,'zscale','log');
set(current_axes, 'yscale','log');
set(current_axes,'xscale','log');

figure(2);
current_axes=gca;
surf(kz_max_graze)
xlabel ('predator vol, \mum^3');ylabel ('prey vol, \mum^3');
zlabel ('kz, umol N/L');
title('kz based on inverse of max ingestion rate');
set(current_axes,'zscale','log');
set(current_axes, 'yscale','log');
set(current_axes,'xscale','log');


figure(3)
current_axes=gca;
surf(kz_h_all)
xlabel ('predator vol, \mum^3');ylabel ('prey vol, \mum^3');
zlabel ('kz, umol N/L');
title('kz based on direct and calculated handling times');
set(current_axes,'zscale','log');
set(current_axes, 'yscale','log');
set(current_axes,'xscale','log');






