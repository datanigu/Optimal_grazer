%%%% Examination of Handling Times %%%%%%%%%%%%%%

%%% Use this code to examine empirical evidence for handling times of
%%% microzooplankton grazing on a variety of prey
%%% Handling times have either been directly measured or are calculated
%%% from the maximum grazing rate

%%% Some of the info in this code came from "Handling_Encounter_Check.m"
%%% from the directory
%%% "Documents/MIT_URI_Postdoc/Grazer_Community_from_First_Principles"
%%% which contains info on both empirical and theoretically-based handling
%%% times

%% Setting the stage
cd('/Users/dtaniguchi/Documents/MIT_URI_Postdoc/Grazer_Community_from_First_Principles/Functional_Groups')
clear all;
close all;
clc;

set(0,'defaultaxesfontsize',30);
set(0,'defaulttextfontsize',30);
set(0,'defaultlinelinewidth',2);

%% Here I'm going to try to fit a curve that is based on handling time
% values 
%**********THESE ARE ACTUAL MEASUREMENTS OF HANDLING TIME, NOT DERIVED FROM MAX GRAZING**********
% Funcitonal form for my paper is:
    %handling time = a * (radius predator/radius prey)^x [maybe should be based on vol?]
    % So I'm going to assuming handling time is some function of the ratio
    % of predator:prey, constrained by the four data points in this study
    % that looked specifically at handling time.  Am only considering
    % handling time of ingested prey
% Values from Boenigk and Arndt (2000, Journal of Eukaryotic Microbiology), for heterotrophic nanoflagellates
% volumes of predators found in table 1, pg. 351, which I will convert to
% radius, length of prey found in table 3, pg. 355
% handling time in table 2, pg. 353, and includes time for contact,
% processing, ingestion, and refractory time
% Handling time is tiem from first contat to resuming of beating of
% flagella [which I think is correct]
% This study includes 1 salt water species and 3 freshwater species
% Temperature is 20 degrees C, pg. 350

predvol_Boenigk = [18.5 25.8 75.4 163];%um^3, Table 1, pg. 351
predr = (3.*predvol_Boenigk./4./pi).^(1/3);% um, Table 3, pg. 355
preyr = [2.6 1.7 1.6 1.4]./2;%um
preyvol_Boenigk = 4/3*pi.*(preyr.^3);%um^3
h_Boenigk = [94.5 33.2 10.1 3.87];%seconds, mean values for ingested particles

t_Boenigk = repmat(20,size(predr));%degrees C

%% Value from Fenchel (1982, MEPS, Ecology of heterotrophic microflagellates: 
% I. Some important forms and their functional morphology, particularly pg. 214
% Predator values and handling time from that source. Prey volume taken from 
% Fenchel (1982, Ecology of heterotrophic micro flagellates: II. Bioenergetics and growth)
predr_Fenchel_h = 7.25/2; %um, for Ochromonas, took midpoint of 7 to 7.5 um, provided diameter on pg. 214
predvol_Fenchel_h = 4/3*pi.*(predr_Fenchel_h.^3);%um^3
h_Fenchel_h = 20;%seconds
preyvol_Fenchel_h = 0.6;%um^3, for prey Pseudomonas (pg. 226)
preyr_Fenchel_h = (3.*preyvol_Fenchel_h./4./pi).^(1/3);% um
prey_C_Fenchel_h = 10e-10; %mg C/cell for Pseudomonas (pg. 226)
t_Fenchel_h = repmat(20,size(h_Fenchel_h));%degrees C, pg. 212, using temperature to which cells transferred every few days rather than stock culture temp

%% Jeong et al., 2005, Aquatic Microbial Ecology 40:133-150
% ***** Actual measurements of handling times (i.e., not based on gmax) for
% mixotrophic dinoflagellates ********
% Examined time it took for mixotrophic dinos to completely engulf an
% unidentified cryptophyte, which had esd of 5.6 um (pg. 139)
% Lingulodinium polyedrum (esd 38.2 um, Table 1) took 73 s
% Prorocentrum donghaiense (esd 13.3, Table 1) took 357 s
% Heterocapsa triquetra (esd 15.0, Table 1) took 549 s
% Prorocentrum micans (esd 26.6, Table 1) took 629 s

% Also examined handling time for L. polyedrum to feed on a variety of prey (pg. 139)
% L. poly feeding on Heterosigma akashiwo (11.5 um esd) was 178 s
% L. poly feeding on Prorocentrum minimum (12.1 um esd) was 293 s
% L. poly feeding on Scrippsiella trochoidea (22.8 um esd) was 756 s
% Above size values taken from Table 1
predrad_Jeong = 0.5.*[13.3 15.0 26.6 38.2 38.2 38.2 38.2];%um
predvol_Jeong = 4/3*pi.*(predrad_Jeong.^3);%um^3
preyrad_Jeong = 0.5.*[5.6 5.6 5.6 5.6 11.5 12.1 22.8];%um
preyvol_Jeong = 4/3*pi.*(preyrad_Jeong.^3);%um^3
h_Jeong = [357 549 629 73 178 293 756];%s
t_Jeong = repmat(20,size(h_Jeong));%degrees C, pg. 138

%% Jeong et al., 2010, Aquatic Microbial Ecology 59: 239-255
% Looked at lots of characteristics of Gymnodinium aureolum, a mixotrophic
% dinoflagellate that feeds via a peduncle (pg. 249)
% To consume a prey cell of the cryptophyte Teleaulax
% sp., beginning with deployment of the peduncle and ending wtih complete
% ingestion, took 256 s (pg. 249).  More than one predator could be attached to one
% prey cell (pg. 249)
% The time it took between deployment of the tow filament [to draw the prey
% closer I think] and to release the peduncle was 27 s (pg. 249)
% When feeding, cell length of predator = 19.8 um and cell width 15.6 um(pg. 246)
%** Does point out that, when fed multiple prey one at a time, growth and
%grazing not significantly correlated with size of prey (Fig. 5)
h_Jeong2010 = 268+27;%s, to account for tow filament and eating through peduncle
t_Jeong2010 = 20;%degrees C, pg. 240
preyrad_Jeong2010 = 0.5*5.6;%um, Table 1
preyvol_Jeong2010 = 4/3*pi*(preyrad_Jeong2010^3);%um^3

predvol_Jeong2010 = 4/3*pi*(19.8/2)*((15.6/2)^2);%um^3
predrad_Jeong2010 = (3/4/pi*predvol_Jeong2010)^(1/3);%um

%% Jeong et al., 2005, Aquatic Microbial Ecology, 38: 249-257
% Looked at feeding of mixotrophic dinoflagellate Gonyaulax polygramma to feed on
% variety of prey, among many other properties.  
% Culturing done at 20 degrees C
% ESD of G. polygramma was 32.5 um (Table 1)
% Fed via both apical horn and sulcus
% Time to feed on Amphidinium carterae (esd 6.6 um, Table 1) via the apical horn was 291 seconds
% and time to feed on same prey via sulcus was 346 s, which were not
% significantly different (pg. 253)
% Mixotroph could only feed on cells <17 um esd, so only fed on sm. things
% (pg. 255)
% This is the first study to show dinoflagellates feed through the apical
% horn and not just the sulcus (pg. 255-256)
% Dino seems to prefer to engult large prey through the sulcus and smaller
% prey through the apical horn (pg. 256)
h_Jeong2005 = (346+291)/2;%s
t_Jeong2005 = 20;%degrees C
preyrad_Jeong2005 = 3.3;%um
preyvol_Jeong2005 = 4/3*pi*(3.3.^3);%um^3
predrad_Jeong2005 = 32.5/2;%um
predvol_Jeong2005 = 4/3*pi*(predrad_Jeong2005.^3);%um^3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% More attempts at emipirically deriving handling times 
% Jeong et al. (2010) for values of heterotrophic dinoflagellates
%Going to use max ingestions rates (ng C/grazer/day), convert to 1/d, and
%then see how they relate to predator:prey size ratios
% Values taken from Table 5
%******** BASED ON MAX GRAZING **************
dinoradius = [6.0 7.8 8.7 13.5 13.5 13.9 15.6 21 21.6 26.8 31 31.8 ...
    36.1 43.5 43.5 45.7 61 73 81]./2; %um

dinomaxingestion = [0.04 2.9 0.8 2.6 4.3 0.8 1.3 2.3 1.3 17.8 19.4 13.6 ...
    11.5 17.1 24.4 17.7 12 5.3 12.1];%ngC/grazer/day
dinopredprey = [0.8 1.3 0.4 2.2 2.2 1.2 1.4 5.3 0.7 1 0.9 2.6 1.4 1.3 1.1 ...
    1.3 1.6 2 2.3];%dimensionless
dinopreyradius = dinoradius./dinopredprey;
dinopreyradius =[8.0 5.9 19.7 6.1 6.1 11.5 11.5 4.0 33.0 26.5 35.2 12.1 ...
    26.5 34.0 37.9 35.2 38.2 26.6 35.2]./2;%um--apparently this are given in the table
dinopreyvol = 4/3*pi.*(dinopreyradius.^3);

dinovol = 4/3*pi.*(dinoradius.^3);
dinoC = (0.76.*(dinovol.^0.819)).*10^-3;%ngC/cell
dinoh = 1./dinomaxingestion.*dinoC*.24.*3600;%s
directengulf = [1 7 12 14 15];%indices just for dinoflagellates that feed by direct engulfment
t_dino = [15 20 20 20 20 20 20 20 20 20 15 20 20 20 20 12 19 19 12];%degrees C

%% Now values from Verity (1991) Limniology and Oceanography, for the ciliates Tintinnopsis dadayi and
% Strobilidium spiralis---NB: I don't think these get used because using
% values from this source but cited in Hansen et al. (1997) below to be
% consistent
% NB: these values are in Hansen et al. (1997) below
cilC = [1.25*10^4 3.2*10^3];%pg C/ciliate, Tintinnopsis and Strobilidium, respectively, pg. 730
Tdadavol = 10.^((log10(cilC(1))-0.639)/0.984);%converting using equation for loricate ciliate given in 
    % Menden-Deuer and Lessard (2000)
Sspirvol = 10.^((log10(cilC(2))-0.168)/0.841);% converting using equation for aloricate ciliate given in 
    % Menden-Deuer and Lessard (2000)
predvol_Verity = [Tdadavol Sspirvol];%um^3
cilradius = (predvol_Verity.*3./4./pi).^(1/3);%um
cilmaxingestion = [0.14  0.14];%1/hour, with data for Tintinnopsis taken from data theif and Fig. 3A
 %For Strobilidum, max ingestions based on single-prey experiments, which didn't differ
 %among different prey species (pg. 734).  The ones measured were for
 % grazing on Isochrysis galbana and Pseudobodo.  For Tintinnopsis, it did vary with prey, so I'm basing it
 %off the prey ingested at the highest rate, Katodinium

cilprey = [8 4.5 3]./2;%um, from Table 2, with the last two for Isochyrsis galbana and Pseudobodo, respecitvely,
    %which both gave similar max ingestion rates for Strobilidium
    
 cilradius = [cilradius(1) cilradius];
 cilh = 1./[0.14 0.14 0.14].*3600;%s
 cilpredprey = cilradius./cilprey;
 t_Verity = 20;%degrees C, pg. 731
 
 %% Values from Hansen et al. (1997) and associated sources listed in table 3 to get more emprircal
 % data to try to figure out handling time
 % Predator volumes and max ingestion rates are listed in Hansen et al. (1997), and size of prey
 % are taken directly from the source listed
%Fenchel (1982) --different from Fenchel (1982) above
 predvol_Fenchel = [75 20 200 190 50 90];%um^3
 preyvol_Fenchel = 0.6.*ones(size(predvol_Fenchel));%um^3, Fenchel
maxingest_Fenchel = [0.86 0.81 0.57 0.8 0.65 0.56];%1/hr
t_Fenchel = repmat(20, size(predvol_Fenchel));%degrees C, Table 3 Hansen et al. 1997


%Eccleston-Parry and Leadbeater, 1994 
% All predators fed same thing
predvol_Eccleston = [54 220 35 75 212 83];%um^3, from table 3 of Hansen et al
preyvol_Eccleston = 0.67.*ones(size(predvol_Eccleston));%um^3, Eccelston-Parry
maxingest_Eccleston = [1.99 0.79 0.69 0.45 0.2 0.3];%1/hr
t_Eccleston = repmat(20, size(predvol_Eccleston));%degrees C, table 3 Hansen et al. 1997

% Values from Andersen, 1988/1989
predvol_Andersen = 40;%um^3
preyvol_Andersen = 0.6;%um^3, Andersen, pg. 525
maxingest_Andersen = 0.33;%1/h
t_Andersen = 15;%degrees C, Hansen et al. 1997

% Andersson, 1989
predvol_Andersson = 50;%um^3
preyvol_Andersson = 4/3*pi*(mean([0.2 0.9])/2)^3;%um^3, Andersson
maxingest_Andersson = 0.42;%1/h
t_Andersson = 20;%degrees C, table 3 Hansen et al. 1997

% Hollen and Boras, 1991
predvol_Hollen = 65.4;
preyvol_Hollen = 0.53;%Hollen, pg. 77
maxingest_Hollen = 0.09;%1/h
t_Hollen = 25;%degrees C, table 3 Hansen et al. 1997

% Strom, 1991
predvol_Strom = 900;%um^3
preyvol_Strom = 4/3*pi*(4.5/2)^3;%um^3, Strom, pg. 104
maxingest_Strom = 0.112;%1/h, looking at Isochrysis galbana, not Synechococcus
t_Strom = 12;%degrees C, table 3 Hansen et al. 1997

% Hansen, 1992
predvol_Hansen92 = 11500;
preyvol_Hansen92 = 2050;%um^3, Table 1 of Hansen 1992
maxingest_Hansen92 = 0.17;%1/hr
t_Hansen92 = 15;%degrees C, table 3 Hansen et al. 1997

% Strom and Buskey, 1993--Already included in data from Jeong et al. (2010)
% Jeong and Latz, 1994--Already included in data from Jeong et al. (2010)

%Jacobson and Anderson, 1993--not good info on prey size, so not going to use
% Heinbokel, 1978--used mixed prey, so not using

% Buskey and Stoecker, 1988
predvol_Buskey = 2.1*10^5;%
maxingest_Buskey = 0.13;%1/h
preyvol_Buskey = 2050;%um^3, taken from HANSEN 1992, not original sources
t_Buskey = 20;%degrees C, table 3 Hansen et al. 1997

% Hansen et al., 1991
predvol_Hansen91 = 9.546*10^4;
maxingest_Hansen91 = 0.215;%1/h
preyvol_Hansen91 = 2050;%Table 2, Hansen et al. 1991
t_Hansen91 = 18;%degrees C, table 3 Hansen et al. 1997


% Jonsson, 1986
predvol_Jonsson = [1.5*10^5 4*10^4];
maxingest_Jonsson = [0.11 0.11];%1/h
preyvol_Jonsson = 4/3*pi.*([9.7 7.9]./2).^3;%Table 3 of Jonsson, 1986
t_Jonsson = [12, 12];%degrees C, table 3 Hansen et al. 1997

% Verity, 1991--one ciliate Strobilidium spiralis had max ingestion when
% consuming both Isochrysis glabana and Pseudobodo, so will use that ciliate twice
% For the other ciliate Tintinnopsis dadayi, i had max ingestion when eating
% Katodinium rotundatum
%NOTE**********These are the same ciliate values as used above, but I'm
%using values from Hansen et al. (1997) and volume ratios to be consistent
%with other sources above, so this coversion is slightly different
predvol_Verity = [2.65*10^4 2.65*10^4 1.13*10^5];
maxingest_Verity = [0.189 0.189 0.155];%1/h
preyvol_Verity = 4/3*pi.*([3 4.5 8]./2).^3;%Table 1 of Verity 1991
t_Verity = repmat(20,size(predvol_Verity));%degrees C, table 3 Hansen et al. 1997

% Bernard and Rassoulzadegan, 1990
% According to my notes in "GrazingRateParameters.m" for the Poulin and
% Franks model parameterization, the grazing rates listed in Hansen et al.
% (1997) are when the ciliate was grazing on Dunaliella minuta, Monochrysis
% lutheri, and Nannochloris sp.  The predator was always Strobidium
% sulcatum.  So I'm going to use the size of those prey items, which are
% given in Table 2 of Bernard and Rassoulzadegan (1990)
predvol_Bernard = [1*10^4 1*10^4 1*10^4];
maxingest_Bernard = [0.481 0.391 0.151];%1/h
preyvol_Bernard = [140 31.12 8.18];% Table 2 Bernard and Rassoulzadegan
t_Bernard = repmat(22,size(predvol_Bernard));

% Verity 1985--didn't measure size of prey, only C and N content, so not
% using

 
%% Values from Hewett (1980) for ciliate grazing on different sized prey
% Estimated handling time from Michaelis-Menten fit to capture rate vs.
% prey density for ciliate Didinium nasutum
% Note that the predator changed size when fed different prey and at
% different prey concentrations, so I'll use the size when fed on the
% particular prey and at high prey concentrations, assuming that is when
% handling time matters (handling time limited rather than encounter
% limited)
% No temperature reported, so no way to temperature correct values and thus
% cannot use them--blah
h_Hewett = [0.916 1.212 1.404];% h, Table 2, pg. 1078, for Paramecium aurelia, P. Jenningsi, and P. multimicronucleatum, respectively
h_Hewett = h_Hewett.*3600;%seconds
predvol_Hewett = [8.21 12.79 10.6].*10^5;%um^3, taken from Table 6, pg. 1080
predr_Hewett = (predvol_Hewett.*3./4./pi).^(1/3);%um

preyvol_Hewett = 4/3*pi.*([85 156 219]./2).^3;%um^3, from pg. 1076
preyr_Hewett = [85 156 219]./2;%um

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Putting together empirical values %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Putting together actual handling time values

% These are all dinos; used below and later in plotting 
h_Jeong_many = [h_Jeong, h_Jeong2010, h_Jeong2005];%s, all dinoflagellates, all direct measurements of handling time
preyr_Jeong_many = [preyrad_Jeong, preyrad_Jeong2010, preyrad_Jeong2005];%um
predr_Jeong_many = [predrad_Jeong, predrad_Jeong2010, predrad_Jeong2005];%um
predvol_Jeong_many = [predvol_Jeong,predvol_Jeong2010, predvol_Jeong2005 ];%um^3
preyvol_Jeong_many = [preyvol_Jeong, preyvol_Jeong2010, preyvol_Jeong2005];%um^3
t_Jeong_many = [t_Jeong, t_Jeong2010, t_Jeong2005];%degrees C

h_vec = [h_Boenigk, h_Fenchel_h, h_Jeong_many];% seconds, handling time
h_preyr = [preyr, preyr_Fenchel_h, preyr_Jeong_many];% um, prey radius
h_predr = [predr, predr_Fenchel_h, predr_Jeong_many];%um, predator radius

h_preyvol = 4/3*pi.*(h_preyr.^3);%um^3, prey volume
h_predvol = 4/3*pi.*(h_predr.^3);%um^3, prey volume

t_handling = [t_Boenigk, t_Fenchel_h, t_Jeong_many];%degrees C


%% Putting together handling times based on max ingestion rate
predvol_Hansen =[predvol_Fenchel predvol_Eccleston predvol_Andersen ...
    predvol_Andersson predvol_Hollen predvol_Strom predvol_Hansen92 predvol_Buskey ...
    predvol_Hansen91 predvol_Jonsson predvol_Verity predvol_Bernard];
predr_Hansen = (predvol_Hansen.*3./4./pi).^(1/3);%um, I assume

preyvol_Hansen = [preyvol_Fenchel preyvol_Eccleston preyvol_Andersen ...
    preyvol_Andersson preyvol_Hollen preyvol_Strom preyvol_Hansen92 preyvol_Buskey ...
    preyvol_Hansen91 preyvol_Jonsson preyvol_Verity preyvol_Bernard];
preyr_Hansen = (preyvol_Hansen.*3./4./pi).^(1/3);%um, I assume

maxingest_Hansen = [maxingest_Fenchel maxingest_Eccleston maxingest_Andersen ...
    maxingest_Andersson maxingest_Hollen maxingest_Strom maxingest_Hansen92 maxingest_Buskey ...
    maxingest_Hansen91 maxingest_Jonsson maxingest_Verity maxingest_Bernard]; 
h_Hansen = 1./maxingest_Hansen.*3600;%s

t_Hansen = [t_Fenchel, t_Eccleston,t_Andersen, t_Andersson, t_Hollen,...
    t_Strom, t_Hansen92, t_Buskey, t_Hansen91, t_Jonsson, t_Verity, t_Bernard]; %degrees C

% h_ingest = [dinoh h_Hansen h_Hewett];%s <-----Hewett didnt' list temperature, so can't use
% predvol_ingest =[dinovol predvol_Hansen predvol_Hewett];%um^3
% preyvol_ingest = [dinopreyvol preyvol_Hansen preyvol_Hewett];%um^3
h_ingest = [dinoh h_Hansen];%s
predvol_ingest =[dinovol predvol_Hansen];%um^3
preyvol_ingest = [dinopreyvol preyvol_Hansen];%um^3

predr_ingest = (3/4/pi.*predvol_ingest).^(1/3);% um
preyr_ingest = (3/4/pi.*preyvol_ingest).^(1/3);%um

t_ingest = [t_dino,t_Hansen];%degrees C

%% Putting together ALL empirical values, based on max ingestion rate AND handling time together
h_all = [h_vec h_ingest];%s

predvol_all = [h_predvol predvol_ingest];
preyvol_all = [h_preyvol preyvol_ingest];

predr_all = (predvol_all.*3./4./pi).^(1/3);%um
preyr_all = (preyvol_all.*3./4./pi).^(1/3);%um

t_all = [t_handling, t_ingest];%degrees C

pred_prey_ratio_all = predr_all./preyr_all;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Temperature correction %%%%%%%%%%%%

% Temperature correction using Metabolic Theory of Ecology
% Based on Chen et al. (2012) L&O and Taniguchi et al. (2014)
% This code is taken from Documents/Poulin and Franks Model/Parameterization/GrazingRateParameters.m 

%Correcting to 20 degrees C = 293.15 K.  Used regression from Fig. 3A. from Chen et al. (2012), slope of 0.67,
%specifically for grazing (frowth has different activation energy value)

boltz = 8.62*10^-5;%eV/Kelvin, Boltzmann's constant
t_corr = 273.15+20;%K, corrected temperature
E = 0.67;%eV, activation energy for max grazing, from Chen et al. (2012)

% Handling time temp corrections
coeff_h = h_vec.*exp(E./boltz./(t_handling+273.15));%coefficient for temperature correction
h_t_corr = coeff_h./exp(E./boltz./t_corr);%seconds

% Max ingestion temp corrections
coeff_ingest = h_ingest.*exp(E./boltz./(t_ingest+273.15));
h_ingest_t_corr = coeff_ingest./exp(E./boltz./t_corr);%s

% All rate temperature corrections
coeff_all = h_all.*exp(E./boltz./(t_all+273.15));
h_all_t_corr = coeff_all./exp(E./boltz./t_corr);%s

% Select groupings for purposes of plotting below
coeff_dino = dinoh.*exp(E./boltz./(t_dino+273.15));
dino_h_t_corr = coeff_dino./exp(E./boltz./t_corr);%s
coeff_Hansen = h_Hansen.*exp(E./boltz./(t_Hansen+273.15));
h_Hansen_t_corr = coeff_Hansen./exp(E./boltz./t_corr);%seconds
coeff_Boenigk=h_Boenigk.*exp(E./boltz./(t_Boenigk+273.15));
h_Boenigk_t_corr = coeff_Boenigk./exp(E./boltz./t_corr);%seconds
coeff_Fenchel_h = h_Fenchel_h.*exp(E./boltz./(t_Fenchel_h+273.15));
h_Fenchel_h_t_corr = coeff_Fenchel_h./exp(E./boltz./t_corr);%s
coeff_Jeong_many = h_Jeong_many.*exp(E./boltz./(t_Jeong_many+273.15));
h_Jeong_many_t_corr = coeff_Jeong_many./exp(E./boltz./t_corr);%s



%% Finding regressions

% Regression for handling time based on direct measurements of handling time
x_direct_hand = [ones(size(h_predr')) log10(h_predr./h_preyr)'];
[param_direct_hand, ci_direct_hand resid_direct_hand outlier_direct_hand ...
    stats_direct_hand]=regress(log10(h_t_corr)',x_direct_hand);
coeff_direct_hand = 10.^param_direct_hand(1);%intercept
exp_direct_hand = param_direct_hand(2);%slope in log space, exponent in linear space
regress_direct_hand = coeff_direct_hand.*(h_predr./h_preyr).^exp_direct_hand;


% Regression for handling time based on inverse of max ingestion rates
x_ingest = [ones(size(predr_ingest')) log10(predr_ingest./preyr_ingest)'];
[param_ingest ci_ingest resid_ingest outlier_ingest stats_ingest]=regress(log10(h_ingest_t_corr)',x_ingest);
coeff_ingest = 10.^ param_ingest(1);
exp_ingest = param_ingest(2);
regress_ingest = coeff_ingest.*(predr_ingest./preyr_ingest).^exp_ingest;
ppratio_ingest = predr_ingest./preyr_ingest;

% Regression for handing time based on all values
x_all = [ones(size(predr_all')) log10(predr_all./preyr_all)'];
[param_all ci_all resid_all outlier_all stats_all] = regress(log10(h_all_t_corr)',x_all);
coeff_all = 10.^param_all(1);
exp_all = param_all(2);
regress_all = coeff_all.*(predr_all./preyr_all).^exp_all;



%% %%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting based on DIRECT MEASUREMENTS OF HANDLING TIME 
% Handling time
figure
loglog(h_predr./h_preyr,h_t_corr,'*','MarkerSize',20)
hold on
loglog(h_predr./h_preyr, regress_direct_hand,'-','LineWidth',10)
legend('direct handling time measurements','linear regression','Location','SouthWest')
legend('boxoff')
xlabel('consumer volume/prey volume')
ylabel('handling time, seconds')


nlx = [ones(size(h_predr)); h_predr./h_preyr]';
lx = [ones(size(h_predr)); log(h_predr./h_preyr)]';
    
[nlparam, nlint, nlr, nlrint, nlstats] = regress(h_t_corr',nlx);
[lparam, lint, lr, lrint, lstats] = regress(log(h_t_corr)', lx);


% Can compare fits for the log data (lparam) and for nonlog data (nlparam)
% and use what has better r^2 value. Below I've used lparam, but not sure
% if that's better
ppratio = predr./preyr;
htestnew = 10^(lparam(1)).*(ppratio).^lparam(2);



% htestnewII = htestnew;
% htestnewII(ppratio<=0.4) = 0;
% disp(['new max handling time for pred:prey size ratios>=0.4 from Boenigk and Arndt 2000 = ',num2str(max(max(htestnewII))/3600/24),' days'])
% disp(['new max handling time for pred:prey size ratios>=0.4 from Boenigk and Arndt 2000= ',num2str(max(max(htestnewII))/3600),' hours'])
% disp(['new max handling time for pred:prey size ratios>=0.4 from Boenigk and Arndt 2000= ',num2str(max(max(htestnewII))),' seconds'])
% disp(['new max handling time for pred:prey size ratios>=1 from Boenigk and Arndt 2000= ',num2str(max(max(htestnew(ppratio>=1)))/3600/24),' days'])
% disp(['new max handling time for pred:prey size ratios>=1 from Boenigk and Arndt 2000= ',num2str(max(max(htestnew(ppratio>=1)))/3600),' hours'])
% disp(['new max handling time for pred:prey size ratios>=1 from Boenigk and Arndt 2000= ',num2str(max(max(htestnew(ppratio>=1)))),' seconds'])
% 
% disp(['new max handling time for pred:prey size ratios>=10 from Boenigk and Arndt 2000= ',num2str(max(max(htestnew(ppratio>=10)))/3600/24),' days'])
% disp(['new max handling time for pred:prey size ratios>=10 from Boenigk and Arndt 2000= ',num2str(max(max(htestnew(ppratio>=10)))/3600),' hours'])
% disp(['new max handling time for pred:prey size ratios>=10 from Boenigk and Arndt 2000= ',num2str(max(max(htestnew(ppratio>=10)))),' seconds'])

%% Plotting handling times JUST based off of MAX INGESTION RATES

figure(11);clf
loglog(predr_ingest./preyr_ingest,h_ingest_t_corr,'+b')
hold on
loglog(predr_ingest./preyr_ingest,regress_ingest,'-','LineWidth',10)
hold on
%Dinos
% loglog(dinovol./dinopreyvol, dinoh,'xr')
loglog(dinoradius./dinopreyradius, dino_h_t_corr,'xr')
hold on
% Mostly ciliates and nanoflagellates
% loglog(predvol_Hansen./preyvol_Hansen, h_Hansen,'db')
loglog(predr_Hansen./preyr_Hansen, h_Hansen_t_corr,'db')

% Ciliate from Hewett
% loglog(predvol_Hewett./preyvol_Hewett, h_Hewett,'+m')
% loglog(predr_Hewett./preyr_Hewett, h_Hewett,'+m')

% Direct engulfment dinoflagellates
directengulf = [1 7 12 14 15];%indices just for dinoflagellates that feed by direct engulfment
% loglog(dinovol(directengulf)./dinopreyvol(directengulf), dinoh(directengulf),'ok','MarkerSize',20)
loglog(dinoradius(directengulf)./dinopreyradius(directengulf), dino_h_t_corr(directengulf),'ok','MarkerSize',20)

hold on
% xlabel('pred vol/preyvol')
xlabel('pred radius/prey radius')
ylabel('handling time, seconds')
legend('all','regression','dinoflagellates','mostly ciliates and flagellates','direct engulf dino','Location','SouthEast')
legend('boxoff')
title('Handling time based off max ingestion')


%% Plotting of all handling times, direct and derived
figure(10);clf
% loglog(pred_prey_ratio_all,h_all,'oc')
loglog(pred_prey_ratio_all,h_all_t_corr,'oc')
xlabel('pred radius/prey radius')
ylabel('handling time, seconds')

hold on
loglog(pred_prey_ratio_all,regress_all,'-k','LineWidth',10)
title('All Handling Times, from inverse max ingestion and direct handling measurements')

%% Comparison of handling time based on max ingestion and actual measurements vs. just max ingestion, for func groups
% Includes delineation among functional groups

% Max ingestion and actual handling time measurements
figure(14);clf
%Flagellates
% loglog([predvol_Boenigk./preyvol_Boenigk predvol_Hansen(1:15)./preyvol_Hansen(1:15), predvol_Fenchel_h./preyvol_Fenchel_h],[h_Boenigk, h_Hansen(1:15),h_Fenchel_h],'+b','MarkerSize',10,'LineWidth',5)
loglog([predr./preyr predr_Hansen(1:15)./preyr_Hansen(1:15), predr_Fenchel_h./preyr_Fenchel_h],...
    [h_Boenigk_t_corr, h_Hansen_t_corr(1:15),h_Fenchel_h_t_corr],'+b','MarkerSize',10,'LineWidth',5)

hold on
%Dinoflagellates
% loglog([dinovol./dinopreyvol predvol_Hansen(16:17)./preyvol_Hansen(16:17),predvol_Jeong_many./preyvol_Jeong_many],[dinoh, h_Hansen(16:17),h_Jeong_many],'xr','MarkerSize',10,'LineWidth',5)
loglog([dinoradius./dinopreyradius predr_Hansen(16:17)./preyr_Hansen(16:17),predr_Jeong_many./preyr_Jeong_many],...
    [dino_h_t_corr, h_Hansen_t_corr(16:17),h_Jeong_many_t_corr],'xr','MarkerSize',10,'LineWidth',5)

hold on
%Ciliates
% loglog([predvol_Hansen(18:end)./preyvol_Hansen(18:end) predvol_Hewett./preyvol_Hewett],[h_Hansen(18:end), h_Hewett],'oc','MarkerSize',10,'LineWidth',5)
loglog([predr_Hansen(18:end)./preyr_Hansen(18:end)],[h_Hansen_t_corr(18:end)],'oc','MarkerSize',10,'LineWidth',5)

% Highlighting actual measurements of handling time
hold on
loglog(h_predr./ h_preyr,h_t_corr,'ok','MarkerSize',20,'LineWidth',3)
% loglog(predvol_Jeong_many./preyvol_Jeong_many,h_Jeong_many,'or','MarkerSize',20)
% loglog(predvol_Boenigk./preyvol_Boenigk,h_Boenigk,'og','MarkerSize',20)
% loglog(predvol_Fenchel_h./preyvol_Fenchel_h, h_Fenchel_h,'ob','MarkerSize',20)
%Regression
hold on
loglog(predr_all./preyr_all,regress_all,'-k','LineWidth',10)
hold on
loglog(h_predr./h_preyr,regress_direct_hand,'-b','LineWidth',5)
hold on
loglog(predr_ingest./preyr_ingest,regress_ingest,'-g','LineWidth',5)

% xlabel('predvol/preyvol')
xlabel('predator radius/prey radius')
ylabel('handling time, seconds')
title('Handling Time')
legend('flagellages','dinoflagellates','ciliates','direct measure of handling','regression all','regression handling','regression inverse max ingest','Location','SouthWest')
legend('boxoff')

% INVERSE OF handling time 
figure(15);clf
%Flagellates
loglog([predr./preyr predr_Hansen(1:15)./preyr_Hansen(1:15), predr_Fenchel_h./preyr_Fenchel_h],...
    [1./h_Boenigk_t_corr, 1./h_Hansen_t_corr(1:15),1./h_Fenchel_h_t_corr],'+b','MarkerSize',10,'LineWidth',5)
hold on
%Dinoflagellates
loglog([dinoradius./dinopreyradius predr_Hansen(16:17)./preyr_Hansen(16:17),predr_Jeong_many./preyr_Jeong_many],...
    [1./dino_h_t_corr, 1./h_Hansen_t_corr(16:17),1./h_Jeong_many_t_corr],'xr','MarkerSize',10,'LineWidth',5)
hold on
%Ciliates
loglog([predr_Hansen(18:end)./preyr_Hansen(18:end)],[1./h_Hansen_t_corr(18:end)],'oc','MarkerSize',10,'LineWidth',5)
% Highlighting actual measurements of handling time
hold on
loglog(h_predr./ h_preyr,1./h_t_corr,'ok','MarkerSize',20,'LineWidth',3)
% xlabel('predvol/prey volume')
xlabel('pred radius/prey radius')
ylabel('1/handling time, 1/seconds')
title('Inverse of handling time')
legend('flagellages','dinoflagellates','ciliates','direct measure of handling','Location','NorthWest')
legend('boxoff')

%%%%%%%%%% Just measure of actual measure of handling time vs pred
%%%%%%%%%% radius:prey radius
figure(22);clf
loglog(h_predr./ h_preyr,h_t_corr,'ok','MarkerSize',20,'LineWidth',3)
hold on
loglog(h_predr./h_preyr,regress_direct_hand,'-b','LineWidth',5)
xlabel('predator radius/prey radius')
ylabel('handling time, seconds')
title('Handling Time')
legend('direct measurements of handling time','linear regression','Location','SouthWest')
legend('boxoff')
xlim([1,10])


%% Color-coordinating values from same source
% NB: STILL BASED ON VOLUME
figure(16);clf
loglog(1./maxingest_Fenchel.*3600,predvol_Fenchel./preyvol_Fenchel,'+k');hold on
loglog(1./maxingest_Eccleston.*3600,predvol_Eccleston./preyvol_Eccleston,'+b');hold on
%loglog(1./maxingest_Andersen.*3600,predvol_Andersen./preyvol_Andersen,'+r');hold on
%loglog(1./maxingest_Andersson.*3600,predvol_Andersson./preyvol_Andersson,'+c');hold on
loglog(1./maxingest_Jonsson.*3600,predvol_Jonsson./preyvol_Jonsson,'+m');hold on
loglog(1./maxingest_Verity.*3600,predvol_Verity./preyvol_Verity,'+y');hold on
loglog(1./maxingest_Bernard.*3600,predvol_Bernard./preyvol_Bernard,'or')
xlabel('predvol/preyvol')
ylabel('handling time, seconds')
title('Handling Time Based on Max Ingestion for Same Sources')

%% Handling time vs. predator size

% All values combined
x_predsize =  [ones(size(predr_all')), log10(predr_all)'];
[param_size, ci_size, resid_size, outlier_size, stats_size] = regress(log10(h_all_t_corr)', x_predsize);
coeff_predsize = 10.^param_size(1);
exp_predsize = param_size(2);
regress_predsize = coeff_predsize.*predr_all.^exp_predsize;

num_h = length(h_vec);
num_ingest = length(h_ingest);

% Regression based on max ingestion
ingest_ind = num_h+1;
x_predsize2 =  [ones(size(predr_ingest')), log10(predr_ingest)'];
[param_size2, ci_size2, resid_size2, outlier_size2, stats_size2] = regress(log10(h_ingest_t_corr)', x_predsize2);
coeff_predsize2 = 10.^param_size2(1);
exp_predsize2 = param_size2(2);
regress_predsize2 = coeff_predsize2.*predr_ingest.^exp_predsize2;

% Regression based on handling time only
x_predsize3 =  [ones(size(h_predr')), log10(h_predr)'];
[param_size3, ci_size3, resid_size3, outlier_size3, stats_size3] = regress(log10(h_t_corr)', x_predsize3);
coeff_predsize3 = 10.^param_size3(1);
exp_predsize3 = param_size3(2);
regress_predsize3 = coeff_predsize3.*h_predr.^exp_predsize3;

% Plotting handling time vs pred size
figure(17);clf
%Flagellates
% loglog([predvol_Boenigk, predvol_Hansen(1:15), predvol_Fenchel_h],[h_Boenigk, h_Hansen(1:15),h_Fenchel_h],'+b','MarkerSize',10,'LineWidth',5)
loglog([predr, predr_Hansen(1:15), predr_Fenchel_h],[h_Boenigk_t_corr, h_Hansen_t_corr(1:15),h_Fenchel_h_t_corr],'+b','MarkerSize',10,'LineWidth',5)

hold on
%Dinoflagellates
loglog([dinoradius, predr_Hansen(16:17),predr_Jeong_many],...
    [dino_h_t_corr, h_Hansen_t_corr(16:17),h_Jeong_many_t_corr],'xr','MarkerSize',10,'LineWidth',5)
hold on
%Ciliates
loglog([predr_Hansen(18:end)],[h_Hansen_t_corr(18:end)],'oc','MarkerSize',10,'LineWidth',5)
% Highlighting actual measurements of handling time
hold on
loglog(h_predr,h_t_corr,'ok','MarkerSize',20,'LineWidth',3)
% Regression for all points
hold on; loglog(predr_all,regress_predsize,'-k','LineWidth',5)
% Regression based on values from max ingestion rate
hold on;loglog(predr_all(ingest_ind:end),regress_predsize2,'-g','LineWidth',5)
% Regression based on values from handling time
hold on;loglog(predr_all(1:num_h),regress_predsize3,'-b','LineWidth',5)

%legend('flagellages','dinoflagellates','ciliates','direct measure of handling','regression','regression based on max ingest','Location','SouthEast')
legend('flagellages','dinoflagellates','ciliates','direct measure of handling','regression all data','regression inverse max ingest','regression handling','Location','SouthEast')

legend('boxoff')

% xlabel('predvol, um^3')
xlabel('pred radius, \mum')
ylabel('handling time, seconds')
title('Handling Time vs. Predator Size')




% Handling time direct measurements, just based on predator, not prey
figure
loglog(h_predr, h_t_corr,'d','MarkerSize',20)
% xlabel('consumer volume, \mum^3')
xlabel('consumer radius, \mum')
ylabel('handling time, seconds')

% Inverse handling time, based on predator only and direct measurements of
% handling
figure
loglog(h_predr, 1./h_t_corr,'d','MarkerSize',20)
% xlabel('consumer volume, \mum^3')
xlabel('consumer radius, \mum')
ylabel('1/handling time, 1/seconds')

% Inverese of handling time vs PRED size (NOT pred:prey ratio), all measurements
figure
%Flagellates
% loglog([predvol_Boenigk, predvol_Hansen(1:15), predvol_Fenchel_h],[1./h_Boenigk, 1./h_Hansen(1:15),1./h_Fenchel_h],'+b','MarkerSize',10,'LineWidth',5)
loglog([predr, predr_Hansen(1:15), predr_Fenchel_h],...
    [1./h_Boenigk_t_corr, 1./h_Hansen_t_corr(1:15),1./h_Fenchel_h_t_corr],'+b','MarkerSize',10,'LineWidth',5)
hold on
%Dinoflagellates
% loglog([dinovol, predvol_Hansen(16:17),predvol_Jeong_many],[1./dinoh, 1./h_Hansen(16:17),1./h_Jeong_many],'xr','MarkerSize',10,'LineWidth',5)
loglog([dinoradius, predr_Hansen(16:17),predr_Jeong_many],...
    [1./dino_h_t_corr, 1./h_Hansen_t_corr(16:17),1./h_Jeong_many_t_corr],'xr','MarkerSize',10,'LineWidth',5)
hold on
%Ciliates
loglog([predr_Hansen(18:end)],[1./h_Hansen_t_corr(18:end)],'oc','MarkerSize',10,'LineWidth',5)
% Highlighting actual measurements of handling time
hold on
loglog(h_predr,1./h_t_corr,'ok','MarkerSize',20,'LineWidth',3)
% xlabel('predvol, \mum^3')
xlabel('pred radius, \mum')
ylabel('1/handling time, 1/seconds')
title('Inverse of handling time')
legend('flagellages','dinoflagellates','ciliates','direct measure of handling','Location','NorthEast')
legend('boxoff')

