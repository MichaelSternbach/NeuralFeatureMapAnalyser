FigureFolder = '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/Figures/';
DataFolder = '/home/michael/Cloud/PhD/MarsupialData/Data/';

%%WallabyH
%%% All Trials
BootstrapSampleFile='TestBig700WallabyH.mat';
JackknifeSampleFile='JackknifeSamplesWallabyH.mat';
PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)

%%% Trial 1-20
BootstrapSampleFile='BootstrapSamples1000WallabyHTrials1-20.mat';
JackknifeSampleFile='JackknifeSamplesWallabyHTrials1-20.mat';
PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)

%%% Trial 1-10
BootstrapSampleFile='BootstrapSamples1000WallabyHTrials1-10.mat';
JackknifeSampleFile='JackknifeSamplesWallabyHTrials1-10.mat';
PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)

% %%WallabyC
% BootstrapSampleFile='TestBig2400WallabyC.mat';
% JackknifeSampleFile='JackknifeSamplesWallabyC.mat';
% PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
% PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)