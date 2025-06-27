clear all
close all
disp('Start OPM analysis')
disp('-------------------------------------------')
%% add path to OPM processing functions 
addpath OPM_Processing/


%% animal parameter
% animal = 'ferret';%'mouse lemur';
% specimen_num = 3;

animal = 'cat_jung';%'mouse lemur';
specimen_num = 5;
resetFilter = true;%true;

%% general parameter data evaluation
getCI = false;
bootstrapsamples = 100;

%% parameter column spacing
smallest_w_mm = 0.1;
w_step_mm = 0.1;
largest_w_mm = 2;
FilterMap = true;

%% parameter filter cut off determination
%splitROI = 4;

lowpass_cutoffs_mm = 0.1:0.01:1.3;
profile_range_mm = 0.1:0.01: 2;
% profile_range_mm = [0.01 1];
% profile_step_mm = 0.01;

%% pw info parameter
do_plotting = false;
%lowpass_cutoffs_mm
beta_rel = 0.5;

%% data folder
AnimalDataFolder = '~/CIDBN/';
ResultDataFolder = ['DataTestGIF' lower(animal) '/' lower(animal) num2str(specimen_num) '/'];
mkdir(ResultDataFolder)

%% plot parameter
FigureFileName = [ResultDataFolder 'ResultOPM_' animal num2str(specimen_num) '.eps'];
border_ColumnSpacing = 20;
position_ColumnSpacing_txt = [10 10];
x_hypercolumns = 1;

%% load data
disp('load data')
[data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,specimen_num,AnimalDataFolder);
data_obj.prepare_samples_array(bootstrapsamples)
disp('-----------------------')
data_obj.activateGIF(true,4)

% %% test data processing
% delete(gcp('nocreate'))
% parpool('local', 2);
% ProcessDataOPM(data,data_info,data_obj.ROI)


%% plot Map
figure;
z1 = data_obj.filter_map(data_obj.read_map()); %map
plot_map(z1,data_obj.ROI,0,1)


% %% test GIF
% data_obj.generateCleanedDataSamplesGIF()
% figure;
% z2 = data_obj.filter_map(data_obj.read_map()); %map
% plot_map(z2,data_obj.ROI,0,1)
% 
% figure;
% plot_mapAbs(abs(z1-z2))
%data_obj = data_handle_corrected(data_info,data,data_obj.ROI);

%% set ROI
data_obj.make_ROI();
ROI = data_obj.ROI;
save([data_path 'exp_info.mat'],'ROI','-append')

%% set filter cut offs
disp('determine filter settings')
splitROI = 2;
% splitROI{1}.x =[138 240];
% splitROI{1}.y =[67 165];
% splitROI{2}.x =[240 345];
% splitROI{2}.y =[67 165];
data_info = getFilterSettings(data_obj,data_info,ResultDataFolder,splitROI,resetFilter,lowpass_cutoffs_mm,profile_range_mm);
disp('-----------------------')
%% calculate column spacing
disp('calculate column spacing')
figure;
[average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,ResultDataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap);
disp('-----------------------')
%% find pinwheels and calc pinwheel density
disp('calculate pinwheel positions and pinwheel density')
PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,ResultDataFolder,newROI,getCI,do_plotting,lowpass_cutoffs_mm,beta_rel);
disp('-----------------------')

%% plot map, pinwheel positions, pinwheel density and column spacing
disp('plot results')

% data_obj.generateCleanedDataSamplesGIF()
f = figure();

z = data_obj.filter_map(data_obj.read_map()); %map
plot_map(z,data_obj.ROI,0,1)

hold on % column spacing
plot([border_ColumnSpacing x_hypercolumns*data_info.pix_per_mm*average_spacing_mm+border_ColumnSpacing],[border_ColumnSpacing border_ColumnSpacing],'-w')
hold on
text(position_ColumnSpacing_txt(1),position_ColumnSpacing_txt(2),[num2str(x_hypercolumns) '\Lambda = ' num2str(round(x_hypercolumns*average_spacing_mm,2)) 'mm'],'Color','white')

hold on %pinwheel positions
plot(PwInfo.PWxList,PwInfo.PWyList,'ow')
hold on;
contour(real(z),[0 0],'white')
hold on;
contour(imag(z),[0 0],'white')

% pinwheel density
title(['pinwheel density=' num2str(round(PwInfo.WeightedPwDensityFixedFilter,2)) '/\Lambda^2'])

% save figure
print(f,'-depsc', FigureFileName)
disp('-----------------------')


% %% get covariance data
% MapCovariances = getMapCovariances(data_obj,ResultDataFolder,'vector',DoFilter,scale);
% NoiseCovariances = getNoiseCovariances(data_obj,ResultDataFolder,'vector',DoFilter,scale);