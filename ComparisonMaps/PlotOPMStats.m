FigureFolder = '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/Figures/';
DataFolder = '/home/michael/Cloud/PhD/MarsupialData/Data/';
addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/AGWolfOPMDataPipeline'

% %%WallabyH
% %%% All Trials
% BootstrapSampleFile='TestBig700WallabyH.mat';
% JackknifeSampleFile='JackknifeSamplesWallabyH.mat';
% PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
% PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)
% 
% %%% Trial 1-20
% BootstrapSampleFile='BootstrapSamples1000WallabyHTrials1-20.mat';
% JackknifeSampleFile='JackknifeSamplesWallabyHTrials1-20.mat';
% PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
% PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)
% 
% %%% Trial 1-10
% BootstrapSampleFile='BootstrapSamples1000WallabyHTrials1-10.mat';
% JackknifeSampleFile='JackknifeSamplesWallabyHTrials1-10.mat';
% PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
% PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)

% %%WallabyC
% %%%All Trials
% BootstrapSampleFile= 'TestBig2400WallabyC.mat';
% JackknifeSampleFile='JackknifeSamplesWallabyC.mat';
% PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
% PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)

% %%% Trial 1-20
% BootstrapSampleFile='BootstrapSamples1000WallabyCTrials1-20.mat';
% JackknifeSampleFile='JackknifeSamplesWallabyCTrials1-20.mat';
% PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
% PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)
% 
% %%% Trial 1-10
% BootstrapSampleFile='BootstrapSamples1000WallabyCTrials1-10.mat';
% JackknifeSampleFile='JackknifeSamplesWallabyHTrials1-10.mat';
% PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
% PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)
% 
% %%Cat BC
% BootstrapSampleFile= 'BootstrapSamples1000CatBCTrials1-10.mat';
% JackknifeSampleFile='JackknifeSamplesCatBCTrials1-10.mat';
% PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
% PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)


% %% Ferret
% experiment_num = 1;%28;
% animal = 'ferret';
% trial_ini=1;
% 
% [data_info,data_path] = info_handle(animal,experiment_num);
% 
% set_blocks = data_info.protocol.blocks;
% trials_to_use = find(set_blocks>0);
% load([data_path,'Processed_2/trial_',num2str(trials_to_use(trial_ini)),'.mat'],'data');
% Bootstrapsamples = load([data_path,'Analyzed_2/characterization/',sprintf('trial_%d_domain_stats',trials_to_use(trial_ini)),'.mat'],'orientation_stats','parameters'); 
% 
% TestPWStatsChepe(animal,experiment_num,data,data_info,data_path,Bootstrapsamples,DataFolder,FigureFolder)



%% Mouse lemur
experiment_num = 1;%28;
animal = 'mouse lemur';
trial_ini=1;


[data_info,data_path] = info_handle(animal,experiment_num);

data_file = ['Processed/',data_info.exp_data.folder,data_info.exp_data.ref];%data_info.ID,'.mat'
load([data_path,data_file],'data');
Bootstrapsamples = 100;
TestPWStatsChepe(animal,experiment_num,data,data_info,data_path,Bootstrapsamples,DataFolder,FigureFolder)
TestOrientationStatsChepe(animal,experiment_num,data,data_info,data_path,Bootstrapsamples,DataFolder,FigureFolder)