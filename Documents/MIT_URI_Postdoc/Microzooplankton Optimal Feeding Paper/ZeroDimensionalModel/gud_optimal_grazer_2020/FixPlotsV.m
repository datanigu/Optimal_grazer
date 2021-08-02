%%%% Fixing Plots created in Fortran for zero-dimensional model

%%% Use this code to plot prettier versions of the microzooplankton
%%% and phytoplankton biomass and the proportion of times each prey size
%%% was grazed, all under different nutrient concentrations

%%% The Fortran code used to create the binary files is "FuncGpsDifGrazeParamV.F"
%%% in the directory
%%% /Users/dtaniguchi/Documents/MIT_URI_Postdoc/Microzooplankton Optimal
%%% Feeding Paper/ZeroDimensionalModel
%%% and the Matlab code that originally plotted the specific plots is FuncGpsPlot.m
%%% (although really only the Fortran code needs to be run)

%% Setting stage
clear all;
close all;
clc;

set(0,'defaultaxesfontsize',30);
set(0,'defaulttextfontsize',30);
set(0,'defaultlinelinewidth',2);

%mainFolder = '/Users/dtaniguchi/Documents/MIT_URI_Postdoc/Microzooplankton Modeling Paper/ZeroDimensionalModel/Ctot_Based_Info_and_Plots';
% mainFolder = '/Users/dtaniguchi/Documents/MIT_URI_Postdoc/Microzooplankton Optimal Feeding Paper/ZeroDimensionalModel/Ctot_Based_Info_and_Plots_II';
mainFolder = '/Users/dtaniguchi/Documents/MIT_URI_Postdoc/Microzooplankton Optimal Feeding Paper/ZeroDimensionalModel/gud_optimal_grazer_2020';

cd(mainFolder)

%% Flags
read_variables_flag = 1;% set to 1 if want to read in variables
all_steps_flag =0;% set to 1 if used all time steps
plot_biomass_flag = 1;%set to 1 if want to plot biomass spectra with different grazers under
    % same nutrient concentrations together
plot_proportion_grazed_flag =1;%plot proportion of times grazed
plot_biomass_allnut_flag = 0;%set to 1 if want to plot together biomass spectra with
    %same grazers under different nutrient concentrations
plot_proportion_allnut_flag = 0;% set to 1 if want to plot, for each nutrient condition,
    %multigrazer systems in which bar plots are grouped side-by-side rather
    %than overlayed as above--only works if plot_biomass_allnut_flag is 1
normalize_flag = 1;% set to 1 if want to plot normalized biomass spectra<----NB This flag normalizes BOTH phyto and microzoopl
scale_flag = 1;%set to 1 if want to have corresponding plots have same x and y lims
biomass_mean_flag = 1;% set to 1 if want mean biomass from last few time points, 0 if just want last time point
save_plots_flag = 0;%set to 1 if want to save plots generated in this code
plot_color_flag = 0;% set to 1 if want plots to be in color, else will be in grayscale
plot_parameters_flag = 0;


%% Nutrient concentrations and grazer sizes used in naming directories
nutconc = [2,1,7];%{'02','06','14'};%which techinically stands for 2.5, 1.5, 0.7 gC/m^3
zsizedir = {'5', '50', '200','Multi'};%which stand for the radius sizes of grazers or for all grazers together (i.e., "Many")
Ctotvec = [2.5 1.5 0.7];%actual total carbon concentrations
colormx = [0 0 1;0.5 0.5 0.5; 0 1 0; 0 0 0];%for colored plots
bwmx = [0.85 0.85 0.85; 0.6 0.6 0.6; 0.3 0.3 0.3; 0 0 0];%colors for black and white (really grayscale)

mmx =['o' 'o' 'o' '+'];%define markers for microzooplankton biomass plots
msizemx = [20 30 40 40];%Marker size when comparing grazers under different nutrient conditions
lmx =[15 25 35 20]; %[15 15 15 8];%line widths for phytoplankton biomass plots
lnutmx = [15 25 35 20];%[10 15 20 20];%line widths for phyto biomass under different nutrient concentrations for each specific grazer
zname = {'05um','50um','200um','zMany'};%grazer systems used to identify directories corresponding to those systems

if plot_color_flag == 1
    plotcolor = colormx;
else
    plotcolor = bwmx;
end

% Plot colors for proportion of times grazed plots
    plotcolorprop(1,:) = plotcolor(1,:);
    plotcolorprop(2,:) = plotcolor(3,:);
    plotcolorprop(3,:) = plotcolor(2,:);%Have to do this because directories are not in monotonic order with respect to grazer size, unlike
        % the proportion of times grazed plot
    plotcolorprop(4,:) =[0 0 0];
% Line style for proportion times grazed for groupings based on dif nut
% conc
lsmx = {':','--','-','-.'};

%% Read in variables
if read_variables_flag == 1
subdirs = dir(mainFolder);

%This filters out all the items in the main folder that are not directories
subdirs(~[subdirs.isdir]) = [];
%And this filters out the parent and current directory '.' and '..'
tf = ismember( {subdirs.name}, {'.', '..'});
subdirs(tf) = [];
numberOfFolders = length(subdirs);
% subdirs_order = {'Ctot_1_5_z05_2018','Ctot_1_5_z50_2018','Ctot_1_5_z200_2018','Ctot_1_5_zmany_200_2018'};
subdirs_order = {'Ctot_1_5_z05_2020','Ctot_1_5_z50_2020','Ctot_1_5_z200_2020','Ctot_1_5_zmany_200_2020'};
 

    for kk = 1:length(subdirs_order)
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

    if all_steps_flag ~=1        
        % tcount = number of time steps actually recorded
        fid=fopen('tcount.bin','r','n');
        tcount= fread(fid,'integer*8');
        fclose(fid);
        tnum = tcount;
        
        % savestep = number of time steps between each saved value
        fid=fopen('savestep.bin','r','n');
        savestep= fread(fid,'integer*8');
        fclose(fid);     

    end %all_steps_flag
    
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
    if all_steps_flag ~=1
        gtimechange=reshape(gtimechangetemp,zmax,tnum_all,fnum);
    else
        gtimechange=reshape(gtimechangetemp,zmax,tnum,fnum);
    end% all_steps_flag
    
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

%% Getting rid of really small values
zchange(zchange<=10^-20) = 10^-20;
pchange(pchange<=10^-20) = 10^-20;

%% Time-related values
start =9e5;%5e5;%7e5;% time point from which to start plotting
endt = tnum;
timestep = 1;%vector for PLOTTING only every other timestep value (e.g., every 10th time step)
% BUT keep in mind that not every time step may have been saved in the
% first place
if biomass_mean_flag ==1 %Determines if want to plot mean biomss or just last time point
    timevec = (start:timestep:endt);
else
    timevec = endt;
end
num_years_plotted = (endt-start)*savestep/3600/24/365*dt;

%% Finding proportion times grazed
gchangemx = NaN(zmax,pnum,fnum);
for j = 1:fnum
    for k = 1:pnum
        for l = 1:zmax
            gchangemx(l,k,j) = length(find(gtimechange(l,start:endt,j)==k));
        end
    end
end

gproportion = gchangemx./(1+endt-start);


%% Renaming matrices and things
%     eval(['p_',subdname,'=pchange;'])
%     eval(['z_',subdname,'=zchange;'])
%     eval(['gprop_',subdname,'=gproportion;'])
    
    % This may be a better, shorter name
    eval(['pchange',num2str(kk),'=pchange;'])
    eval(['zchange',num2str(kk),'=zchange;'])
    eval(['gprop',num2str(kk),'=gproportion;'])
    eval(['zr',num2str(kk),'=zrmx(find(zrmx~=0));'])
    eval(['pr',num2str(kk),'=pr;'])
    eval(['Ctot',num2str(kk),'=Ctot;'])
    
    clear zchange pchange gtimechange gchangemx gproportion

    
    end%kk, loop through all directories  
end%read_variables_flag


%% Plotting steady state solutions in 2-D plots
if plot_biomass_flag == 1
    figcounter =0;%for the biomass plots under each nutrient condition
    
%     for ll = 1:numel(nutconc)
ll = 2;
        figcounter = figcounter+1;
        ccount = 0;% for defining color of plots
        
%         % Find subdirectories with a specific nutrient concentration
%         ind1 = strfind({subdirs.name}, num2str(nutconc(ll)));
%         indtemp =find(~cellfun(@isempty,ind1));
%         
%         % Sorting indtemp by z values so they are in the order of zname
%         for iji = 1:length(indtemp)
%             %eval(['ztemp(',num2str(iji),') = zr',num2str(ind(iji)),';'])
%             ind2 = strfind({subdirs.name}, char(zname(iji)));
%             indind = find(~cellfun(@isempty,ind2));
%             ind(iji) = intersect(indtemp,indind);
%         end%iji

        for jj = 1:4%going through each data set with the same nutrient concentration
            ccount = ccount+1;
            %Renaming values so can more easily plot
            eval(['zchange = zchange',num2str(jj),';'])
            eval(['pchange = pchange',num2str(jj),';'])
            eval(['gprop = gprop',num2str(jj),';'])
            eval(['zr = zr',num2str(jj),';'])
            eval(['pr = pr',num2str(jj),';'])
            
%%     Microzooplankton biomass spectrum plot--NB: what the mean is of is dealt with earlier in defining "timevec"
            figure(figcounter+10)
            set(gcf,'Position',[43,1,1238,705])
            if normalize_flag == 1 %Want to normalize biomass spectrm
              loglog(2.*zr.*10^6, mean(zchange(1:length(zr),timevec),2)./(zr.*zsizewidth),mmx(ccount),'LineWidth',6,'MarkerSize',30,'Color',plotcolor(ccount,:));hold on                                
                ylabel('grazer biomass, gC m^-^2')
                if scale_flag==1 % Want all plots to be scale similarly
                    %ylim([0.2*10^-16 5*10^5])
                    ylim([0 5*10^6])
                end               
            else %DO NOT want to normalize spectrum
                loglog(2.*zr.*10^6, mean(zchange(1:length(zr),timevec),2),mmx(ccount),'LineWidth',6,'MarkerSize',30,'Color',plotcolor(ccount,:));hold on
                ylabel('grazer biomass, gC m^-^3')
                if scale_flag==1
                    ylim([0.2*10^-20 5])
                end
            end%if statement for normalized biomass    
            xlabel('grazer diameter, \mum')   
            eval(['title(''Microzooplankton Biomass for Nut Conc ',num2str(Ctotvec(ll)),' g C m^-^3'')'])
            
%% Phytoplankton biomass spectrum--NB: what the mean is of is dealt with earlier in defining "timevec"
            figure(figcounter);
            set(gcf,'Position',[43,1,1238,705])
            %set(gcf,'Position',[1,1,1002,574])
            if normalize_flag == 1 %Want to normalize biomass spectrm
                loglog(2.*pr(1:end-1).*10^6, mean(pchange(1:end-1,timevec),2)./diff(pr),'-','LineWidth',lmx(ccount),'MarkerSize',30,'Color',plotcolor(ccount,:));hold on                
                ylabel('phytoplankton biomass, gC m^-^2')
                if scale_flag==1 % Want all plots to be scale similarly
                    %ylim([10^-20 10^6])
                    ylim([0 5*10^6])
                end               
            else %DO NOT want to normalize spectrum
                semilogy(2.*pr.*10^6, mean(pchange(:,timevec),2),'-','LineWidth',lnmx(ccount),'MarkerSize',30,'Color',plotcolor(ccount,:));hold on
                ylabel('phytoplankton biomass, gC m^-^3')
                if scale_flag==1
                    ylim([10^-14 10^8])
                end
            end%if statement for normalized biomass            
            xlabel('phytoplankton diameter, \mum')  
            %eval(['title(''Phytoplankton Biomass for Dissolved Conc ',num2str(Ctotvec(ll)),' g C m^-^3'')'])
            eval(['title(''Phytoplankton Biomass for Dissolved Conc ',num2str(Ctotvec(ll)),' g C m^-^3'')'])

        end%jj, which goes through each data set with same nut conc

        % Saving phytoplankton plots
        figure(figcounter)
           legend(zsizedir)
            legend('boxoff')
        if (save_plots_flag == 1 && scale_flag == 1)
            eval(['savefig(''Phytoplankton_Biomass_Rescaled_',num2str(Ctotvec(ll)),'Ctot.fig'')'])
            eval(['print(gcf,''Phytoplankton_Biomass_Rescaled_',num2str(Ctotvec(ll)),'Ctot.eps'',''-depsc2'')'])
        elseif (save_plots_flag == 1 && scale_flag ~=1)
            eval(['savefig(''Phytoplankton_Biomass_',num2str(Ctotvec(ll)),'Ctot.fig'')'])
            eval(['print(gcf,''Phytoplankton_Biomass_',num2str(Ctotvec(ll)),'Ctot.eps'',''-depsc2'')'])
        end
        
        %Saving microzooplankton plots
        figure(figcounter+10)
             legend(zsizedir,'Location','NorthWest')
            legend('boxoff')
            if (save_plots_flag == 1 && scale_flag == 1)
                eval(['savefig(''Grazer_Biomass_Rescaled_',num2str(Ctotvec(ll)),'Ctot.fig'')'])
                eval(['print(gcf,''Grazer_Biomass_Rescaled_',num2str(Ctotvec(ll)),'Ctot.eps'',''-depsc2'')'])
            elseif (save_plots_flag == 1 && scale_flag ~=1)
                eval(['savefig(''Grazer_Biomass_',num2str(Ctotvec(ll)),'Ctot.fig'')'])
                eval(['print(gcf,''Grazer_Biomass_',num2str(Ctotvec(ll)),'Ctot.eps'',''-depsc2'')'])
            end

        
%     end%ll, going through each nutrient concentration
end%plot_biomass_flag
    
%% Plotting proporiton of times each prey size class was grazed
if plot_proportion_grazed_flag == 1
    prop_fig_counter = 20;%for the figure numbers
    ii=2;
%     for ii = 1:numel(nutconc)
%         % Find subdirectories with a specific nutrient concentration
%         ind2 = strfind({subdirs.name}, num2str(nutconc(ii)));
%         ind3temp =find(~cellfun(@isempty,ind2));
%         
% % Sorting indtemp by z values so they are in the order of zname
%         for iji = 1:length(ind3temp)
%             ind2temp = strfind({subdirs.name}, char(zname(iji)));
%             indind = find(~cellfun(@isempty,ind2temp));
%             ind3(iji) = intersect(ind3temp,indind);
%         end%iji        
%         
        ccount = 0;
        for kk = 1:4%ind3
            ccount = ccount+1;
            %Renaming values so can more easily plot
            eval(['gprop = gprop',num2str(kk),';'])
            eval(['zr = zr',num2str(kk),';'])
            eval(['pr = pr',num2str(kk),';'])
            prop_fig_counter = prop_fig_counter+1;
            figure(prop_fig_counter);clf
                set(gcf,'Position',[43,1,1238,705])
                 numgrazers = length(find(sum(gprop,2)~=0)); %number of grazers
                if numgrazers==1 %meaning there is only one grazer
                bar((1:length(pr)),gprop(1,:),'FaceColor',plotcolor(ccount,:))
                eval(['title(''Proportion Times Grazed for Grazer ',num2str(2.*zr.*10^6),' \mum Nut Conc ',num2str(Ctotvec(ii)),' g C m^-^3'')'])
                else
                  bartemp = bar((1:length(pr)),gprop(1:numgrazers,:)',1,'grouped');
                  for i = 1:numgrazers
                      %bar((1:length(pr)),gprop(i,:),1/i,'FaceColor',plotcolorprop(i,:));hold on
%                       set(bartemp(i),'FaceColor',plotcolorprop(i,:));
                        set(bartemp(i),'FaceColor',plotcolor(i,:));
                      eval(['title(''Proportion Times Grazed for Multiple Grazers under Nut Conc ',num2str(Ctotvec(ii)),' g C m^-^3'')'])
                      legendtext{i} = num2str(2*zr(i)*10^6);
                  end%i
                  legend(legendtext)
                  legend('boxoff')
                end%what to plot based on number of grazers
                xlabel('prey diameter, \mum')
                
                set(gca,'XTick',1:pnum)
                set(gca,'XTickLabel','')
                set(gca,'XTickLabel',2.*pr*10^6)
                ylabel('proprotion times grazed')
                ylim([0 1])
                xlim([0 pnum+1])
                if save_plots_flag == 1%This is always scaled
                    if numgrazers == 1
                    eval(['savefig(''Proportion_Times_Grazed_Z',num2str(2.*zr*10^6),'_',num2str(Ctotvec(ii)),'Ctot_Rescaled.fig'')'])
                    eval(['print(gcf,''Proportion_Times_Grazed_Z',num2str(2.*zr*10^6),'_',num2str(Ctotvec(ii)),'Ctot_Rescaled.eps'',''-depsc2'')'])
                    else
                    eval(['savefig(''Proportion_Times_Grazed_Zmany_',num2str(Ctotvec(ii)),'Ctot_Rescaled.fig'')'])
                    eval(['print(gcf,''Proportion_Times_Grazed_Zmany_',num2str(Ctotvec(ii)),'Ctot_Rescaled.eps'',''-depsc2'')'])

                    end%if numgrazers
                end    
            
            
        end%kk, going through each directory with a specific nutrient concentration
%     end%ii, going through each nutrient concentration
end%plot_proportion_flag

%% Examining biomass values
p1_tot = sum([sum(mean(pchange1(:,timevec),2)),sum(mean(zchange1(1:length(zr),timevec),2))]);
p2_tot = sum([sum(mean(pchange2(:,timevec),2)),sum(mean(zchange2(1:length(zr),timevec),2))]);
p3_tot = sum([sum(mean(pchange3(:,timevec),2)),sum(mean(zchange3(1:length(zr),timevec),2))]);
p4_tot = sum([sum(mean(pchange4(:,timevec),2)),sum(mean(zchange4(1:length(zr),timevec),2))]);


%% Plotting same grazers, different nutrient concentrations together on same plot (for the most part)
if plot_biomass_allnut_flag == 1
fcounter = 0;%for biomass plots for each grazer
lcounter = 0;%for the legend

for gg = 1:numel(zname) %going through each grazer system
    fcounter = fcounter+1;
    ccount = 0;%for defining color of plots
    pcount1 = 0;%for arrays for proportion of times grazed
    pcount2 = 0;
    % Find subdirectories with a specific grazer size/name
        indI = strfind({subdirs.name}, char(zname(gg)));
        indII =find(~cellfun(@isempty,indI));
    for ggg = indII%going through each nutrient concentration
        ccount = ccount+1;
        lcounter = lcounter+1;
            %Renaming values so can more easily plot
            eval(['zchange = zchange',num2str(ggg),';'])
            eval(['pchange = pchange',num2str(ggg),';'])
            eval(['gprop = gprop',num2str(ggg),';'])
            eval(['zr = zr',num2str(ggg),';'])
            eval(['pr = pr',num2str(ggg),';'])
            eval(['Ctot = Ctot',num2str(ggg),';'])
            
%%     Microzooplankton biomass spectrum plot--NB: what the mean is of is dealt with earlier in defining "timevec"
numgrazers = length(find(sum(gprop,2)~=0)); %number of grazers
if numgrazers==1 %meaning there is only one grazer
    figure(51)
    set(gcf,'Position',[43,1,1238,705])
    if normalize_flag == 1 %Want to normalize biomass spectrm
        semilogy(Ctot, mean(zchange(1:length(zr),timevec),2)./(zr.*zsizewidth),mmx(ccount),'LineWidth',3,'MarkerSize',msizemx(gg),'Color',plotcolor(gg,:));hold on
        ylabel('grazer biomass, g C m^-^2')
        if scale_flag==1 % Want all plots to be scale similarly
            %ylim([0.2*10^-16 5*10^5])
            ylim([0 5*10^2])
        end
    else %DO NOT want to normalize spectrum
        semilogy(Ctot, mean(zchange(1:length(zr),timevec),2),mmx(ccount),'LineWidth',5,'MarkerSize',msizemx(gg),'Color',plotcolor(gg,:));hold on
        ylabel('grazer biomass, gC m^-^3')
        if scale_flag==1
            ylim([0.2*10^-20 5])
        end
    end%if statement for normalized biomass
    xlabel('total nutrient concentration, g C m^-^3')
    legend1{lcounter} = [num2str(2*zr*10^6),' ',num2str(Ctot)];
    xlim([0.5 3])
else %more than one grazer
    figure(52);
    set(gcf,'Position',[43,1,1238,705])
    for iii = 1:numgrazers
        if normalize_flag == 1 %Want to normalize biomass spectrm
            semilogy(Ctot, mean(zchange(iii,timevec),2)./(zr(iii).*zsizewidth),'+','LineWidth',iii*4,'MarkerSize',msizemx(iii),'Color',...%plotcolorprop(iii,:));hold on
                plotcolor(iii,:));hold on
            ylabel('grazer biomass, g C m^-^2')
            if scale_flag==1 % Want all plots to be scale similarly
                %ylim([0.2*10^-16 5*10^5])
                ylim([0 5*10^6])
            end
        else %DO NOT want to normalize spectrum
            semilogy(Ctot, mean(zchange(iii,timevec),2),'+','LineWidth',iii*4,'MarkerSize',msizemx(iii),'Color',...%plotcolorprop(iii,:));hold on
                plotcolor(iii,:));hold on
            ylabel('grazer biomass, gC m^-^3')
            if scale_flag==1
                ylim([0.2*10^-20 5])
            end
        end%if statement for normalized biomass
        legend2{iii} = num2str(2*zr(iii)*10^6);
    end%ii, going through each grazer size
    xlabel('total nutrient concentration, g C m^-^3')   
    xlim([0.5 3])
end%if numgrazers==1

%% Phytoplankton biomass spectrum--NB: what the mean is of is dealt with earlier in defining "timevec"
            figure(fcounter+40);
            set(gcf,'Position',[43,1,1238,705])
            if normalize_flag == 1 %Want to normalize biomass spectrm
                loglog(2.*pr(1:end-1).*10^6, mean(pchange(1:end-1,timevec),2)./diff(pr),'-',...
                    'LineWidth',lnutmx(ccount),'LineStyle',lsmx{ccount},'Color',...%plotcolorprop(gg,:));hold on
                    plotcolor(gg,:));hold on
                ylabel('phytoplankton biomass, gC m^-^2')
                if scale_flag==1 % Want all plots to be scale similarly
                    %ylim([10^-20 10^6])
                    ylim([0 5*10^6])
                end               
            else %DO NOT want to normalize spectrum
                semilogy(2.*pr.*10^6, mean(pchange(:,timevec),2),'-',...
                    'LineWidth',lnutmx(ccount),'LineStyle',lsmx{ccount},'Color',...%plotcolorprop(gg,:));hold on
                    plotcolor(gg,:));hold on
                ylabel('phytoplankton biomass, gC m^-^3')
                if scale_flag==1
                    ylim([10^-14 10^8])
                end
            end%if statement for normalized biomass            
            xlabel('phytoplankton diameter, \mum')  
            legend3{ccount} = [num2str(Ctot),' gC m-3'];

            %% Setting up arrays and matrices for eventual plotting of proportion of times grazed
            if plot_proportion_allnut_flag == 1
                if numgrazers==1 %meaning there is only one grazer
                    pcount1 = pcount1+1;
                    gpropN(pcount1,:,gg) = gprop(1,:);
                    CtotN(pcount1) = Ctot;
                else
                    pcount2 = pcount2+1;
                    for ki = 1:numgrazers
                        gpropNMulti(pcount2,:,ki) =gprop(ki,:);
                        CtotNMulti(pcount2) = Ctot;
                    end%ki
                end
            end            
            
    end%ggg, going through each nutrient concentration  
    
%     %% Renaming proportion times grazed matrices and arrays
%     if plot_proportion_allnut_flag == 1
%     eval(['gpropN',zname(gg),'=gpropN;'])
%     eval(['gpropNMulti',zname(gg),'=gpropNMulti;'])
%     end
    
    %% Saving figures
    % Saving phytoplankton plots
    figure(fcounter+40)
    legend(legend3)
    legend('boxoff')
    eval(['title(''Phytoplankton Biomass for Grazer ',zname{gg},''')'])
    if (save_plots_flag == 1 && scale_flag == 1)
        eval(['savefig(''Phytoplankton_Biomass_Rescaled_',zname{gg},'_System.fig'')'])
        eval(['print(gcf,''Phytoplankton_Biomass_Rescaled_',zname{gg},'_System.eps'',''-depsc2'')'])
    elseif (save_plots_flag == 1 && scale_flag ~=1)
        eval(['savefig(''Phytoplankton_Biomass_',zname{gg},'_System.fig'')'])
        eval(['print(gcf,''Phytoplankton_Biomass_',zname{gg},'_System.eps'',''-depsc2'')'])
    end
    
end%gg, going through each grazer system

%% Title and saving microzooplankton plots
        %Saving microzooplankton plots
        figure(51)
        legend(legend1,'Location','NorthEast')
        legend('boxoff')
        title('Microzooplankton Biomass vs. Total Nutrients, Single Grazer')
        if (save_plots_flag == 1 && scale_flag == 1)
            savefig('Grazer_Biomass_Rescaled_vs_Ctot_Single.fig')
            print(gcf,'Grazer_Biomass_Rescaled_vs_Ctot_Single.eps','-depsc2')
        elseif (save_plots_flag == 1 && scale_flag ~=1)
            savefig('Grazer_Biomass_Not_Scaled_vs_Ctot_Single.fig')
            print(gcf,'Grazer_Biomass_Not_Scaled_vs_Ctot_Single.eps','-depsc2')
        end
        
        figure(52)
        legend(legend2,'Location','NorthEast')
        legend('boxoff')
        title('Microzooplankton Biomass vs. Total Nutrients, Multi Grazer')
        if (save_plots_flag == 1 && scale_flag == 1)
            savefig('Grazer_Biomass_Rescaled_vs_Ctot_Multi.fig')
            print(gcf,'Grazer_Biomass_Rescaled_vs_Ctot_Multi.eps','-depsc2')
        elseif (save_plots_flag == 1 && scale_flag ~=1)
            savefig('Grazer_Biomass_Not_Scaled_vs_Ctot_Multi.fig')
            print(gcf,'Grazer_Biomass_Not_Scaled_vs_Ctot_Multi.eps','-depsc2')
        end
        
        %% Proportion of times each prey size class was grazed
        if plot_proportion_allnut_flag == 1
            fpcounter = 0;
            for ij = 1:3%going through each grazer system (Only go through three because the single and mutli grazer systems are kept separate)
                                
                % Redefinig a few things
                fpcounter = fpcounter+1;
%                 eval(['gpropN = gpropN',zname(gg),';'])
%                 eval(['gpropNMutli = gpropNMulti',zname(gg),';'])
                
                %Plotting single grazer systems
                figure(fpcounter+70)
                set(gcf,'Position',[43,1,1238,705])
                bar1 = bar((1:length(pr)),gpropN(:,:,ij)','grouped');
                for i = 1:length(nutconc)
                    %set(bar1(i),'FaceColor',plotcolor(ij,:),'LineStyle',lsmx{i},'LineWidth',lnutmx(i));
                    %set(bar1(i),'FaceColor',plotcolorprop(ij,:),'LineStyle',lsmx{i},'LineWidth',i);
                    set(bar1(i),'FaceColor',plotcolor(ij,:),'LineStyle',lsmx{i},'LineWidth',i)
                    legendtext1{i} = num2str(CtotN(i));
                end%i
                eval(['title(''Proportion Times Grazed for Single Grazer ',zname{ij},' Under Dif Nut Regimes'')'])
                legend(legendtext1)
                legend('boxoff')
                xlabel('prey diameter, \mum')
                set(gca,'XTick',1:pnum)
                set(gca,'XTickLabel','')
                set(gca,'XTickLabel',2.*pr*10^6)
                ylabel('proprotion times grazed')
                ylim([0 1])
                xlim([0 pnum+1])
                if save_plots_flag == 1%This is always scaled
                    eval(['savefig(''Proportion_Times_Grazed_',zname{ij},'_Dif_Nut_Rescaled.fig'')'])
                    eval(['print(gcf,''Proportion_Times_Grazed_',zname{ij},'_Dif_Nut_Rescaled.eps'',''-depsc2'')'])
                end%save plots flag
                
                % Plotting multigrazer systems
                figure(fpcounter+80)
                set(gcf,'Position',[43,1,1238,705])
                bar2 = bar((1:length(pr)),gpropNMulti(:,:,ij)','grouped');
                 for i = 1:length(nutconc)
                    %set(bar2(i),'FaceColor',plotcolorprop(ij,:),'LineStyle',lsmx{i},'LineWidth',i);
                    set(bar2(i),'FaceColor',plotcolor(ij,:),'LineStyle',lsmx{i},'LineWidth',i);
                    legendtext2{i} = num2str(CtotNMulti(i));
                end%i
                eval(['title(''Proportion Times Grazed for Grazer ',zname{ij},' In Multi Grazer System'')'])
                legend(legendtext2)
                legend('boxoff')
                xlabel('prey diameter, \mum')
                set(gca,'XTick',1:pnum)
                set(gca,'XTickLabel','')
                set(gca,'XTickLabel',2.*pr*10^6)
                ylabel('proprotion times grazed')
                ylim([0 1])
                xlim([0 pnum+1])
                if save_plots_flag == 1%This is always scaled
                    eval(['savefig(''Proportion_Times_Grazed_',zname{ij},'_Dif_Nut_Multi_Rescaled.fig'')'])
                    eval(['print(gcf,''Proportion_Times_Grazed_',zname{ij},'_Dif_Nut_Multi_Rescaled.eps'',''-depsc2'')'])
                end%save plots flag                            
                
            end%ij, for each grazer system
        end%plot_proportion_allnut_flag
            
end%plot_biomass_allnut_flag


%% Variation of parameters with size
if plot_parameters_flag ==1

%zr =[20 40 80 ]*10^-6;%1.0e-03 *[0.0010    0.0020    0.0040    0.0080    0.0160    0.0320    0.0640    0.1280    0.2560];
zr=[1:0.5:500].*10^-6;
v = 20.*zr;%2.^[0:length(zr)-1];
%pr=  1.0e-03 *[0.0004  0.0025 0.005 0.01 0.02 0.0400    0.0800    0.1600    0.3200    0.6400];
pr = [0.5:0.5:256].*10^-6;



% Phytoplankton
mu0 = 1.36/24/3600;% 1/s, value taken from my manuscript
muexp = -0.16;% dimensionless, it's the exponent
mu = mu0.*(2*pr*10^6).^muexp;


ks0 = 0.33*1000*106/16/10^6*12;%0.33*106/16*12;%
ksexp=0.48;
ks = ks0.*(2*pr*10^6).^ksexp;
ksII = (0.33*(2*pr*10^6).^ksexp)*106/16*12/1000;

%% Grazer parameters
addpath('/Users/dtaniguchi/Documents/MIT_URI_Postdoc/Grazer_Community_from_First_Principles')
zv = 4/3*pi*(zr.^3);%m^3, grazer volume
zC_cell = 0.216*((zv*10^18).^0.939);% pg C/cell, using the most general equation from Menden-Deuer and Lessard, 2000
zC_cell = zC_cell.*10^-12;%gC/cell
pv = 4/3*pi*(pr.^3);%m^3, prey volume
pC_cell = 0.216*((pv*10^18).^0.939);% pg C/cell, using the equation from Menden-Deuer and Lessard, 2000
pC_cell = pC_cell.*10^-12;%gC/cell
pC_cell_mx = repmat(pC_cell,length(zr),1);
zC_cell_mx = repmat(zC_cell',1,length(pr));
gamma = 0.7;

% Max grazing
gmax = 33.96.*(( (2.*zr) .*10^6).^-0.66)./24/3600;%1/s, values taken from my parameterization manuscript
gmax_mx = repmat(gmax',1,length(pr));

%Handling time
h = handlingII(gmax_mx, pC_cell_mx, zC_cell_mx,zr,pr,3);


% Capture probability
cexp = 3;%exponent for capture probability
capture_mx = capture_efficiency(pr, zr, cexp);

% Encounter kernel
sal = 0.00035; %kg/kg = 35 psu
tmp = 15;%C, temperature

etaW = 4.2844*10^(-5.0)+(0.157*(tmp+64.993)^2.0-91.296)^(-1.0);
eta_D = etaW*...
    ( 1+ (1.541+(1.998*10^(-2.0))*tmp-(9.52*10^(-5.0))*tmp^2.0)* ...
    (sal) + ...
    ((7.974-7.561*10^(-2.0))*tmp+...
    (4.724*10^(-4.0))*tmp^2.0)*(sal^2.0))*1000;
             

%Diffusion
Dz = ((tmp+273.15)*(1.38*10^-20))./(6*pi*eta_D.*zr);
Dp = ((tmp+273.15)*(1.38*10^-20))./(6*pi*eta_D.*pr);

for i = 1:length(pr)
    for j = 1:length(zr)
 betad(j,i) = 4*pi*(Dp(i) + Dz(j)).*(pr(i)+zr(j));
 beta5(j,i) = pi*v(j)*(zr(j) + pr(i))^2;%
    end
end

beta = betad + beta5;

% Respiration
eta_R = eta_D*10^-3;% viscosity for respiration from swimming
efficiency = 0.01;%efficiency for ciliary propulsion

z0R = 10.^( 0.75*log10(zv*10^18) - 4.09 );%nL O2/cell/hour, basal respiration rate--need because using assimilation efficiency
R0 = z0R./(zv.*10^18) * (1/10^9) * (1/0.08206/(tmp + 273.15)) *12*(10^12) / 0.216/3600;%1/sec, basal respiration rate
    Rv = 3*pi.*zr.*2.*(v.^2).*eta_R/20.2*(10^-3)/(tmp+273.15)/0.082*12./zC_cell/efficiency;
R = R0 + Rv;%1/s

% Mortality
zm0 = 0.05/24/3600;% 
zm = zm0 + zm0.*v;

%Total loss
loss_tot = R + zm;


%% Plotting parameters

zrind = NaN(3,1);
zrind(1) = find(zr == 2.5*10^-6);
zrind(2) = find(zr == 25*10^-6);
zrind(3) = find(zr == 250*10^-6);
legendname = {'5 \mum consumer','50 \mum consumer','500 \mum consumer'};
prxlim = [0.75 105];
zrxlim = [0.75 500];

%% Plotting Phytoplankton Parameters
figure;
  set(gcf,'Position',[43,1,1238,705])
loglog(pr.*10^6,mu,'-k','LineWidth',20)

xlabel('prey radius, \mum')
ylabel('max phytopl growth rate, 1/s')
title('Phytoplankton growth rate')
xlim(prxlim)
%     eval(['savefig(''Phyto Growth.fig'')'])
%     eval(['print(gcf,''Phyto Growth.eps'',''-depsc2'');'])


figure
  set(gcf,'Position',[43,1,1238,705])
loglog(pr.*10^6,ksII,'-k','LineWidth',20)
xlabel('prey radius, \mum')
ylabel('phyto half sat, gC/s')
title('Phytoplankton half sat for nutrients')
xlim(prxlim)
%     eval(['savefig(''Phyto Half Sat.fig'')'])
%     eval(['print(gcf,''Phyto Half Sat.eps'',''-depsc2'');'])
%     
    nutrients = [0.1:0.2:100];
    pr1 = find(pr == 4e-6);
    pr30 = find(pr == 8e-6);
    mumax1 = mu(pr1).*(nutrients./(nutrients + ks(pr1)));
    mumax30 = mu(pr30).*(nutrients./(nutrients + ks(pr30)));
    
    figure
     set(gcf,'Position',[43,1,1238,705])
     plot(nutrients,mumax1,'-k','LineWidth',7)
     hold on
     plot(nutrients, mumax30,'-k','LineWidth',15)
     xlabel('nutrients')
     ylabel('phyto growth rate, 1/s')
     title('Phyto Growth for 4 and 8 \mum Radius Prey')
%         savefig('Phyto Growth for 4 and 8 um Phyto.fig')
 

%% Plotting Grazer parameters
figure
  set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,gmax','-','Color',[0.5 0.5 0.5],'LineWidth',7)
xlabel('grazer radius, \mum')
ylabel('max grazing rate, 1/s')
xlim(zrxlim)
%     eval(['savefig(''Grazer Max Grazing.fig'')'])
%     eval(['print(gcf,''Grazer Max Grazing.eps'',''-depsc2'');'])


figure
  set(gcf,'Position',[43,1,1238,705])
%loglog(zr.*10^6,R,'-','Color',[0.5 0.5 0.5],'LineWidth',7)
loglog(zr.*10^6,R,'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('respiration, 1/s')
xlim(zrxlim)
   eval(['savefig(''Grazer Respiration.fig'')'])
    eval(['print(gcf,''Grazer Respiration.eps'',''-depsc2'');'])
  
    figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,Rv,'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('swimming respiration, 1/s')
xlim(zrxlim)

figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,R0,'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('basal respiration, 1/s')
xlim(zrxlim)
    
figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,zr.*40.*10^6,'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('swimming speed, \mum/s')
xlim(zrxlim)

figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,zm0+zm0*v*1000,'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('mortality, 1/s')
xlim(zrxlim)

figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,zm0*v,'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('swimming mortality, 1/s')
xlim(zrxlim)

figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,0.08291*(zr.*10^6).^2.817.*10^-12./12*1000,'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('cell content, mmol C/cell')
xlim(zrxlim)

figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,loss_tot,'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('total loss, 1/s')
xlim(zrxlim)

figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,1./h(:,1).*0.1./(0.1+zC_cell_mx(:,1)./h(:,1)./capture_mx(:,1)./beta(:,1)),'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('total intake, \mum/s')
xlim(zrxlim)

%% Plots that depend on both grazer and prey
figure
  set(gcf,'Position',[43,1,1238,705])
loglog(pr.*10^6, h(zrind(1),:),'-b','LineWidth',7);hold on
loglog(pr.*10^6, h(zrind(2),:),'-g','LineWidth',7); hold on
loglog(pr.*10^6, h(zrind(3),:),'-r','LineWidth',7)
legend(legendname,'Location','NorthWest')
legend('boxoff')
xlabel('prey radius, \mum')
ylabel('handling time, s')
xlim(prxlim)
%    eval(['savefig(''Handling Time.fig'')'])
%     eval(['print(gcf,''Handling Time.eps'',''-depsc2'');'])
    
    figure % Handling time for just one predator
  set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6, h(:,10),'-','Color',[0.5 0.5 0.5],'LineWidth',7);hold on
xlabel('grazer radius, \mum')
ylabel('handling time, s')
xlim(zrxlim)
%    eval(['savefig(''Handling Time One Predator.fig'')'])
   
     figure % Handling time for one predator, several prey
  set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6, h(:,2),'-k','LineWidth',10);hold on%1 um diameter prey
loglog(zr.*10^6,h(:,97),'-','Color',[0.5 0.5 0.5],'LineWidth',10);hold on%20 um diameter prey
loglog(zr.*10^6,h(:,497),'-','Color',[0.8 0.8 0.8],'LineWidth',10);hold on%100 um diamter prey
xlabel('grazer radius, \mum')
ylabel('handling time, s')
xlim(zrxlim)  
legend('1 \mum prey','20 \mum prey','100 \mum prey')
legend('boxoff')


%Max grazing rate, i.e. 1/h
figure % Handling time for one predator, several prey
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6, 1./h(:,2),'-k','LineWidth',10);hold on%1 um diameter prey
loglog(zr.*10^6,1./h(:,97),'-','Color',[0.5 0.5 0.5],'LineWidth',10);hold on%20 um diameter prey
loglog(zr.*10^6,1./h(:,497),'-','Color',[0.8 0.8 0.8],'LineWidth',10);hold on%100 um diamter prey
xlabel('grazer radius, \mum')
ylabel('max grazing rate, 1/s')
xlim(zrxlim)  
legend('1 \mum prey','20 \mum prey','100 \mum prey')
legend('boxoff')
   

figure% capture probability
  set(gcf,'Position',[43,1,1238,705])
loglog(pr.*10^6,capture_mx(zrind(1),:),'-b','LineWidth',10)
hold on
loglog(pr.*10^6,capture_mx(zrind(2),:),'-g','LineWidth',10)
hold on
loglog(pr.*10^6,capture_mx(zrind(3),:),'-r','LineWidth',10)
legend(legendname,'Location','SouthWest')
legend('boxoff')
xlabel('prey radius, \mum')
ylabel('capture probability')
xlim(zrxlim)
%    eval(['savefig(''Capture Probability.fig'')'])
%     eval(['print(gcf,''Capture Probability.eps'',''-depsc2'');'])
    
    figure %capture probability for one predator
  set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,capture_mx(:,10),'-k','Color',[0.5 0.5 0.5],'LineWidth',10)
xlabel('grazer radius, \mum')
ylabel('capture probability')
xlim(zrxlim)
%    eval(['savefig(''Capture Probability One Predator.fig'')'])

   figure %capture probability for one PREY
  set(gcf,'Position',[43,1,1238,705])
loglog(pr.*10^6,capture_mx(10,:),'-k','LineWidth',10)
xlabel('prey radius, \mum')
ylabel('capture probability')
xlim(prxlim)
%    eval(['savefig(''Capture Probability One Prey.fig'')'])
   
   figure%Capture probability for one predator, three prey
   set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6, capture_mx(:,2),'-k','LineWidth',10);hold on%1 um diameter prey
loglog(zr.*10^6,capture_mx(:,97),'-','Color',[0.5 0.5 0.5],'LineWidth',10);hold on%20 um diameter prey
loglog(zr.*10^6,capture_mx(:,497),'-','Color',[0.8 0.8 0.8],'LineWidth',10);hold on%100 um diamter prey
xlabel('grazer radius, \mum')
ylabel('capture probability')
xlim(zrxlim)
legend('1 \mum prey','20 \mum prey','100 \mum prey','Location','SouthEast')
legend('boxoff')

   
figure%Encounter kernel
  set(gcf,'Position',[43,1,1238,705])
loglog(pr.*10^6,beta(zrind(1),:),'-b','LineWidth',10)
hold on
loglog(pr.*10^6,beta(zrind(2),:),'-g','LineWidth',10)
hold on
loglog(pr.*10^6,beta(zrind(3),:),'-r','LineWidth',10)
legend(legendname,'Location','NorthWest')
legend('boxoff')
xlabel('prey radius, \mum')
ylabel('encounter kernel, m^3/s')
xlim(prxlim)
%    eval(['savefig(''Encounter Kernel.fig'')'])
%     eval(['print(gcf,''Encounter Kernel.eps'',''-depsc2'');'])
%     
    figure %Encounter Kernel one predator
  set(gcf,'Position',[43,1,1238,705])
loglog(pr.*10^6,beta(zrind(1),:),'-','Color',[0.5 0.5 0.5],'LineWidth',10)
xlabel('prey radius, \mum')
ylabel('encounter kernel, m^3/s')
xlim(prxlim)
%    eval(['savefig(''Encounter Kernel One Predator.fig'')'])
   
    
    figure %Encounter Kernel  different prey
  set(gcf,'Position',[43,1,1238,705])
  loglog(zr.*10^6,beta(:,2),'-k','LineWidth',10);hold on
loglog(zr.*10^6,beta(:,97),'-','Color',[0.5 0.5 0.5],'LineWidth',10);hold on
loglog(zr.*10^6,beta(:,497),'-','Color',[0.8 0.8 0.8],'LineWidth',10);hold on
xlabel('grazer radius, \mum')
ylabel('encounter kernel, m^3/s')
xlim(zrxlim)
   legend('1 \mum prey','20 \mum prey','100 \mum prey','Location','SouthEast')
legend('boxoff')
    

% kz = Qz/c/h/beta
figure
set(gcf,'Position',[43,1,1238,705])
  loglog(zr.*10^6,10^-9/12.*0.08291*(zr.*10^6)'.^2.817./h(:,2)./capture_mx(:,2)./beta(:,2),'-k','LineWidth',10);hold on
loglog(zr.*10^6,10^-9/12.*0.08291*(zr.*10^6)'.^2.817./h(:,97)./capture_mx(:,97)./beta(:,97),'-','Color',[0.5 0.5 0.5],'LineWidth',10);hold on
loglog(zr.*10^6,10^-9/12.*0.08291*(zr.*10^6)'.^2.817./h(:,497)./capture_mx(:,497)./beta(:,497),'-','Color',[0.8 0.8 0.8],'LineWidth',10);hold on
xlabel('grazer radius, \mum')
ylabel('grazing half sat, mmol C/m^3')
xlim(zrxlim)
   legend('1 \mum prey','20 \mum prey','100 \mum prey','Location','NorthWest')
legend('boxoff')


figure%kz with no capture prob
set(gcf,'Position',[43,1,1238,705])
  loglog(zr.*10^6,0.08291*zr'.^2.817./h(:,2)./beta(:,2),'-k','LineWidth',10);hold on
loglog(zr.*10^6,0.08291*zr'.^2.817./h(:,97)./beta(:,97),'-','Color',[0.5 0.5 0.5],'LineWidth',10);hold on
loglog(zr.*10^6,0.08291*zr'.^2.817./h(:,497)./beta(:,497),'-','Color',[0.8 0.8 0.8],'LineWidth',10);hold on
xlabel('grazer radius, \mum')
ylabel('grazing half sat, m^3/s')
xlim(zrxlim)
   legend('1 \mum prey','20 \mum prey','100 \mum prey','Location','SouthEast')
legend('boxoff')

% Grazing rate 
% kz = Qz/c/h/beta
figure
%Assume pbiomass = 0.1;
pbiomass = 0.1;
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,0.7./h(:,2)*pbiomass/(pbiomass+0.08291*zr'.^2.817./h(:,2)./capture_mx(:,2)./beta(:,2)),'-k','LineWidth',10);hold on
loglog(zr.*10^6,0.7./h(:,97)*pbiomass/(pbiomass+08291*zr'.^2.817./h(:,97)./capture_mx(:,97)./beta(:,97)),'-','Color',[0.5 0.5 0.5],'LineWidth',10);hold on
loglog(zr.*10^6,0.7./h(:,497)*pbiomass/(pbiomass+0.08291*zr'.^2.817./h(:,497)./capture_mx(:,497)./beta(:,497)),'-','Color',[0.8 0.8 0.8],'LineWidth',10);hold on
xlabel('grazer radius, \mum')
ylabel('grazing rate (not max), m^3/s')
xlim(zrxlim)
   legend('1 \mum prey','20 \mum prey','100 \mum prey','Location','SouthEast')
legend('boxoff')

% Total intake
figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,1./h(:,1).*0.1./(0.1+zC_cell_mx(:,1)./h(:,1)./capture_mx(:,1)./beta(:,1)),'-k','LineWidth',20); hold on
loglog(zr.*10^6,1./h(:,97).*0.1./(0.1+zC_cell_mx(:,97)./h(:,97)./capture_mx(:,97)./beta(:,97)),'-','LineWidth',20,'Color',[0.5 0.5 0.5]);hold on
loglog(zr.*10^6,1./h(:,497).*0.1./(0.1+zC_cell_mx(:,497)./h(:,497)./capture_mx(:,497)./beta(:,497)),'-','LineWidth',20);
xlabel('grazer radius, \mum')
ylabel('total intake, \mum/s')
xlim(zrxlim)

% Total intake for several different prey sizes
figure
set(gcf,'Position',[43,1,1238,705])
loglog(zr.*10^6,1./h(:,1).*0.1./(0.1+zC_cell_mx(:,1)./h(:,1)./capture_mx(:,1)./beta(:,1)),'-k','LineWidth',20)
xlabel('grazer radius, \mum')
ylabel('total intake, \mum/s')
xlim(zrxlim)


   %% Handling time in 3 d (time vs. predator size vs. prey size)
   phytoex = [1:10];
   grazex = [1:20];
   pC_cell = 0.216*(phytoex*4/3*pi).^2.8;
   zC_cell = 0.216*(grazex*4/3*pi).^2.8;
   [pmx, zmx] = meshgrid(phytoex, grazex);
   handling_ex = handlingII(gmax, pC_cell, zC_cell,grazex,phytoex,3);
   handling_ex2 = (3600)*(zmx./(pmx.^5)).^(-1.272);
   
   figure
   subplot(2,1,1)
   surf(phytoex,grazex,log(handling_ex./3600));colorbar
   zlabel('handling')
   ylabel('predator')
   xlabel('prey')
   subplot(2,1,2)
   surf(phytoex,grazex,log(1./handling_ex.*3600));colorbar
   zlabel('graze max')
   ylabel('predator')
   xlabel('prey')
   
      figure;
   subplot(2,1,1)
   surf(phytoex,grazex,log(handling_ex2./3600));colorbar
   zlabel('handling_sq')
   ylabel('predator')
   xlabel('prey')
   subplot(2,1,2)
   surf(phytoex,grazex,log(1./handling_ex2.*3600));colorbar
   zlabel('graze max_sq')
   ylabel('predator')
   xlabel('prey')

end %plot_parameters_flag

