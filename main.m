clear all
close all

%% add path to OPM processing functions 
addpath OPM_Processing/


%% animal parameter
animal = 'ferret';
specimen_num = 6;

%% general parameter data evaluation
getCI = false;
bootstrapsamples = 100;

%% parameter column spacing
smallest_w_mm = 0.1;
w_step_mm = 0.05;
largest_w_mm = 1.5;
FilterMap = true;

%% parameter filter cut off determination
resetFilter = false;
lowpass_cutoffs_mm = 0.1:0.01:0.8;
profile_range_mm = [0.1 2];
profile_step_mm = 0.1;

%% pw info parameter
do_plotting = false;
%lowpass_cutoffs_mm
beta_rel = 0.5;

%% data folder
AnimalDataFolder = '~/CIDBN1/';
ResultDataFolder = ['Data/' lower(animal) '/' lower(animal) num2str(specimen_num) '/'];
mkdir(ResultDataFolder)

%% plot parameter
FigureFileName = ['ResultOPM_' animal num2str(specimen_num) '.eps'];
border_ColumnSpacing = 20;
position_ColumnSpacing_txt = [10 10];
x_hypercolumns = 1;

%% load data
[data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,specimen_num,AnimalDataFolder);
data_obj.prepare_samples_array(bootstrapsamples)

%% set filter cut offs
data_info = getFilterSettings(data_obj,data_info,ResultDataFolder,resetFilter,lowpass_cutoffs_mm,profile_range_mm,profile_step_mm);

%% calculate column spacing
[average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,ResultDataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap);

%% find pinwheels and calc pinwheel density
PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,ResultDataFolder,newROI,getCI,do_plotting,lowpass_cutoffs_mm,beta_rel);


%% plot map, pinwheel positions, pinwheel density and column spacing

f = figure();

z = data_obj.filter_map(data_obj.read_map()); %map
plot_map(z)

hold on % column spacing
plot([border_ColumnSpacing x_hypercolumns*data_info.pix_per_mm*average_spacing_mm+border_ColumnSpacing],[border_ColumnSpacing border_ColumnSpacing],'-w')
hold on
text(position_ColumnSpacing_txt(1),position_ColumnSpacing_txt(2),[num2str(x_hypercolumns) '\Lambda = ' num2str(round(x_hypercolumns*average_spacing_mm,2)) 'mm'],'Color','white')

hold on %pinwheel positions
plot(PwInfo.PWxList,PwInfo.PWyList,'ow')

% pinwheel density
title(['pinwheel density=' num2str(round(PwInfo.WeightedPwDensityFixedFilter,2)) '/\Lambda^2'])

% save figure
print(f,'-depsc', FigureFileName)


% %% get covariance data
% Covariances = getMapCovariances(data_obj,ResultDataFolder,'vector',DoFilter,scale);