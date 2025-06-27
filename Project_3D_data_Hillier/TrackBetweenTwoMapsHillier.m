% clear all
close all
disp('Start OPM analysis')
disp('-------------------------------------------')
%% add path to OPM processing functions 
addpath OPM_Processing/


%% animal parameter
animal = 'Cat_PwTracking';


%% general parameter data evaluation
getCI = false;
bootstrapsamples = 100;


%% data folder
AnimalDataFolder = '/home/michael/Cloud/PhD/data/DataHillier/PwTracking/';
ResultDataFolder = ['Data/' lower(animal) '/'];
mkdir(ResultDataFolder)

%% plot parameter
FigureFileName = [ResultDataFolder 'ResultOPM_' animal '.eps'];
border_ColumnSpacing = 20;
position_ColumnSpacing_txt = [10 10];
x_hypercolumns = 1;

%% load data1
disp('load data1')
FileName1 = '135329_phase_maps_params_masks.mat';
FileData1 = load([AnimalDataFolder FileName1]);
disp('-----------------------')


%% prepare data_info
data_info.animal = animal;
data_info.ID = 'Test1';
data_info.stim_order = FileData1.orientations_in_order;%[0,22.500000000000000,45,67.500000000000000,90,1.125000000000000e+02,135,1.575000000000000e+02];
data_info.pix_per_mm = 16.666666666666670;
data_info.pixels_per_mm = data_info.pix_per_mm;


%% load data2
disp('load data2')
FileName2 = '141752_phase_maps_params_masks.mat';
FileData2 = load([AnimalDataFolder FileName2]);
disp('-----------------------')

%% prepare data
data1 = makeDataFromPhaseMap(FileData1.per_cycle_phase_maps,data_info.stim_order);
data2 = makeDataFromPhaseMap(FileData2.per_cycle_phase_maps,data_info.stim_order);

%% prepare ROI
ROI = FileData1.mid_mask.*FileData1.mid_mask;
data_info.field_size_pix = size(ROI);%[444,408];


%% Make Data Objs
 
data_obj1 = data_handle_corrected(data_info,data1,ROI);
data_obj1.prepare_samples_array(bootstrapsamples)

data_obj2 = data_handle_corrected(data_info,data2,ROI);
data_obj2.prepare_samples_array(bootstrapsamples)



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
DataFile = [ResultDataFolder 'PwData_' animal '.mat'];
save(DataFile,'data_obj1','data_obj2','data_info','StatsOPM1','StatsOPM2','TrackingResults')

%% plot maps
disp('plot maps')

FigureFileName1 = [ResultDataFolder 'Pw1_' animal '.eps'];
SizesCI1 = plotPinwheelCI(z1,data_info,data_obj1.ROI,StatsOPM1.pinwheel_stats,TrackingResults.ini(:,1),FigureFileName1);

FigureFileName2 = [ResultDataFolder 'Pw2_' animal '.eps'];
SizesCI2 = plotPinwheelCI(z2,data_info,data_obj2.ROI,StatsOPM2.pinwheel_stats,TrackingResults.ini(:,2),FigureFileName2);


FigureFileName12 = [ResultDataFolder 'PwMovement_' animal '.eps'];
plotPinwheelMovement(z1,data_info,ROI,StatsOPM1.pinwheel_stats,StatsOPM2.pinwheel_stats,TrackingResults,FigureFileName12)

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

FigureFileName4 = [ResultDataFolder 'PwMovementStatsComparison_' animal '.eps'];
print(f4,'-depsc', FigureFileName4)



function StatsOPM = getStatsOPM(data_obj)
    
    %% get OPM CIs
    alpha = 0.05;
    [StatsOPM.CI_angle,StatsOPM.CI_Abs,StatsOPM.ROI] = getCI(data_obj,alpha,'bca');

    %% get PW stats
    tracker_obj = pinwheel_tracker;
    simple_track=true;
    [StatsOPM.pinwheel_stats,StatsOPM.pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
end

