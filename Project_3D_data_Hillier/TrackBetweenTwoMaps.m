% clear all
close all
disp('Start OPM analysis')
disp('-------------------------------------------')
%% add path to OPM processing functions 
addpath OPM_Processing/


%% animal parameter
animal = 'mouse lemur';
SpecimenNumber = 3;

%% general parameter data evaluation
getCI = false;
bootstrapsamples = 100;


%% data folder
AnimalDataFolder = '~/CIDBN/';
ResultDataFolder = ['Data/' lower(animal) '/' lower(animal) num2str(SpecimenNumber) '/'];
mkdir(ResultDataFolder)

%% plot parameter
FigureFileName = [ResultDataFolder 'ResultOPM_' animal num2str(SpecimenNumber) '.eps'];
border_ColumnSpacing = 20;
position_ColumnSpacing_txt = [10 10];
x_hypercolumns = 1;

%% load data1
disp('load data1')
[data_obj1,data1,data_info1,BloodVesselImg1]=loadMouseLemurData(SpecimenNumber,AnimalDataFolder,1,1);
data_obj1.prepare_samples_array(bootstrapsamples)
disp('-----------------------')

%% load data2
disp('load data2')
[data_obj2,data2,data_info2,BloodVesselImg2]=loadMouseLemurData(SpecimenNumber,AnimalDataFolder,2,1);
data_obj2.prepare_samples_array(bootstrapsamples)
disp('-----------------------')



%% analyse data1
StatsOPM1 = getStatsOPM(data_obj1);

%% analyse data2
StatsOPM2 = getStatsOPM(data_obj2);

%% track between two mean maps
tracker_obj = pinwheel_tracker;
z1 = data_obj1.filter_map(data_obj1.read_map());
z2 = data_obj2.filter_map(data_obj2.read_map());
TrackingResults = tracker_obj.interpolate(z1,z2,data_obj1.ROI);

%% Save data
DataFile = [ResultDataFolder 'PwData_' animal num2str(SpecimenNumber) '.mat'];
save(DataFile,'data_obj1','data_obj2','data_info1','data_info2','StatsOPM1','StatsOPM2','TrackingResults')
%load(DataFile,'data_obj1','data_obj2','data_info1','data_info2','StatsOPM1','StatsOPM2','TrackingResults')

%% plot maps
disp('plot maps')

FigureFileName1 = [ResultDataFolder 'Pw1_' animal num2str(SpecimenNumber) '.eps'];
SizesCI1 = plotPinwheelCI(z1,data_info1,data_obj1.ROI,StatsOPM1.pinwheel_stats,TrackingResults.ini(:,1),FigureFileName1);

FigureFileName2 = [ResultDataFolder 'Pw2_' animal num2str(SpecimenNumber) '.eps'];
SizesCI2 = plotPinwheelCI(z2,data_info1,data_obj2.ROI,StatsOPM2.pinwheel_stats,TrackingResults.ini(:,2),FigureFileName2);

%% compare pw tracking stats
PwTrackingStats = comparePwTrackingStats(StatsOPM1.pinwheel_stats,StatsOPM2.pinwheel_stats,SizesCI1,SizesCI2,TrackingResults);

%% plot Pw Movement Stats

f3 = figure();
t = tiledlayout(1,2);
nexttile;
bins = linspace(-180,180,10);
h_direction=histogram(PwTrackingStats.directions,bins);
xlabel('directions [°]')
ylabel('Absolute Frequency')
xlim([-180 180])
xticks(bins)

nexttile;
histogram(PwTrackingStats.distances)
xlabel('distances [px]')
ylabel('Absolute Frequency')

FigureFileName3 = [ResultDataFolder 'PwMovementStats_' animal '.eps'];
print(f3,'-depsc', FigureFileName3)


f4 = figure();
t = tiledlayout(2,2);
nexttile;
plot(1-PwTrackingStats.prob1,PwTrackingStats.distances,'*')
xlabel('1-PW Prob. OPM1')
ylabel('distance Pw movement [px]')

nexttile;
plot(1-PwTrackingStats.prob2,PwTrackingStats.distances,'*')
xlabel('1-PW Prob. OPM2')
ylabel('distance Pw movement [px]')

nexttile;
plot(PwTrackingStats.SizeCI1,PwTrackingStats.distances,'*')
xlabel('PW CI size OPM1 [px]')
ylabel('distance Pw movement [px]')

nexttile;
plot(PwTrackingStats.SizeCI2,PwTrackingStats.distances,'*')
xlabel('PW CI size OPM2 [px]')
ylabel('distance Pw movement [px]')

FigureFileName4 = [ResultDataFolder 'PwMovementStatsComparison_' animal num2str(SpecimenNumber) '.eps'];
print(f4,'-depsc', FigureFileName4)


function [data_obj,data,data_info,BloodVesselImg]=loadMouseLemurData(SpecimenNumber,AnimalDataFolder,day,experiment)
    
    %% input
    animal = 'mouse lemur';
    if nargin <3
        day = 1;
    end
    if nargin<4
        experiment = 1;
    end
    
    %% get data info
    [data_info,data_path] = info_handle(animal,SpecimenNumber,AnimalDataFolder);
    if isfield(data_info,'pix_per_mm')
        data_info.pixels_per_mm = data_info.pix_per_mm;
    else
        data_info.pix_per_mm = data_info.pixels_per_mm;
    end

    %% load data
    DataFile = [data_path 'Processed/' data_info.exp_data(day).folder data_info.exp_data(day).binocular{experiment}];
    load(DataFile,'data')
    
    %% make blodvessel image
    BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order.binocular);
    
    %% make data object
    data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
    data_obj.apply_LSM()
end

function StatsOPM = getStatsOPM(data_obj)
    
    %% get OPM CIs
    alpha = 0.05;
    [StatsOPM.CI_angle,StatsOPM.CI_Abs,StatsOPM.ROI] = getCI(data_obj,alpha,'bca');

    %% get PW stats
    tracker_obj = pinwheel_tracker;
    simple_track=true;
    [StatsOPM.pinwheel_stats,StatsOPM.pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
end

