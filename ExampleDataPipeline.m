clear all
close all
disp('Start OPM analysis')
disp('-------------------------------------------')
%% add path to OPM processing functions 
addpath OPM_Processing/


%% animal parameter
animal = 'cat';
specimen_num = 1;

%% load data
disp('load data')
load('~/Cloud/PhD/data/DataHillier/124002_phase_maps_params_masks.mat')
% [data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,specimen_num,AnimalDataFolder);
% data = ??? %% [pix_x,pix_y,orientations, trials]
ROI = top_mask;%ones(data_info.field_size_pix);

%% data Parameter
data_info.animal = animal;
data_info.ID = 'Test1';
data_info.field_size_pix = size(ROI);%[444,408];
data_info.stim_order = orientations_in_order;%[0,22.500000000000000,45,67.500000000000000,90,1.125000000000000e+02,135,1.575000000000000e+02];
data_info.pix_per_mm = 16.666666666666670;
data_info.pixels_per_mm = data_info.pix_per_mm;

%% prepare data
data = makeDataFromPhaseMap(per_cycle_phase_maps,data_info.stim_order);


%% input data
bootstrapsamples = 100;
data_obj = data_handle_corrected(data_info,data,ROI);
data_obj.prepare_samples_array(bootstrapsamples)
disp('-----------------------')



%% general parameter data evaluation
CI = false;
alpha = 0.05;
%% parameter column spacing
smallest_w_mm = 0.1;
w_step_mm = 0.05;
largest_w_mm = 1.5;
FilterMap = true;

%% parameter filter cut off determination
resetFilter = false;
lowpass_cutoffs_mm = 0.1:0.01:0.8;
profile_range_mm = [0.01 1.5];
profile_step_mm = 0.01;
% profile_range_mm = [0.01 1];
% profile_step_mm = 0.01;

%% pw info parameter
do_plotting = false;
%lowpass_cutoffs_mm
beta_rel = 0.5;

%% data folder
ResultDataFolder = ['Data/' lower(animal) '/' lower(animal) num2str(specimen_num) '/'];
mkdir(ResultDataFolder)

%% plot parameter
FigureFileName = [ResultDataFolder 'ResultOPM_' animal num2str(specimen_num) '.eps'];
border_ColumnSpacing = 20;
position_ColumnSpacing_txt = [10 10];
x_hypercolumns = 1;

%% contour plots
grey = [0.6 0.6 0.6];
color_contur = grey;%'grey';
linewidth = .1;
set(gca,'Fontsize',20)

[YROI,XROI] = find(data_obj.ROI);
[Xmin, Xmax] = findBorders(XROI);
[Ymin, Ymax] = findBorders(YROI);


%% determine filter cut offs
disp('determine filter settings')
data_info = getFilterSettings(data_obj,data_info,ResultDataFolder,resetFilter,lowpass_cutoffs_mm,profile_range_mm,profile_step_mm);
disp('-----------------------')
%% calculate column spacing
disp('calculate column spacing')
[average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,ResultDataFolder,smallest_w_mm,largest_w_mm,w_step_mm,CI,FilterMap);
disp('-----------------------')
%% find pinwheels and calc pinwheel density
disp('calculate pinwheel positions and pinwheel density')
PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,ResultDataFolder,newROI,CI,do_plotting,lowpass_cutoffs_mm,beta_rel);
disp('-----------------------')

%% plot map, pinwheel positions, pinwheel density and column spacing
disp('plot results')
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

%% plot pw PDF

sigma = 0.05;
PWxList = PwInfo.pinwheel_stats.x(~isnan(PwInfo.pinwheel_stats.x));
PWyList = PwInfo.pinwheel_stats.y(~isnan(PwInfo.pinwheel_stats.y));

local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PWxList, PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
local_pw_dens = local_pw_dens./sum(local_pw_dens,'all').*sum(PwInfo.pinwheel_stats.probability,'all');

% ax = nexttile;
figure;
plot_mapAbs(local_pw_dens,'Pinwheel Prob. Density',max(local_pw_dens(data_obj.ROI),[],'all'),min(local_pw_dens(data_obj.ROI),[],'all'),data_obj.ROI)

hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
hold on
contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
hold on
contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
hold on
for ii = 1:size(PwInfo.pinwheel_stats.x,1)
    text(PwInfo.pinwheel_stats.x(ii,1)+2,PwInfo.pinwheel_stats.y(ii,1),num2str(PwInfo.pinwheel_stats.probability(ii)),'Color',color_contur,'FontSize',5)
end

xlim([Xmin Xmax])
ylim([Ymin Ymax])

yticks([])
xticks([])
print('-depsc', [ResultDataFolder 'CI_PW.eps'])
close

%% get CI OPM
[CI_angle,CI_Abs,ROI] = getCI(data_obj,alpha,'bca');

figure
plot_mapAbs(CI_Abs,' uncertainty selectivity [<selectivity>]',max(CI_Abs,[],'all'),min(CI_Abs,[],'all'),ROI)
print('-depsc', [ResultDataFolder 'CI_Abs.eps'])
close

figure
plot_mapAbs(CI_angle,'uncertainty orientation [°]',180,0,ROI)
print('-depsc', [ResultDataFolder 'CI_Angle.eps'])
close


% %% get covariance data
% MapCovariances = getMapCovariances(data_obj,ResultDataFolder,'vector',DoFilter,scale);
% NoiseCovariances = getNoiseCovariances(data_obj,ResultDataFolder,'vector',DoFilter,scale);