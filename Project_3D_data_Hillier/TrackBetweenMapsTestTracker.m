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

%% get data info
[data_info,data_path] = info_handle(animal,SpecimenNumber,AnimalDataFolder);
if isfield(data_info,'pix_per_mm')
    data_info.pixels_per_mm = data_info.pix_per_mm;
else
    data_info.pix_per_mm = data_info.pixels_per_mm;
end

N=length(data_info.exp_data);

% for ii=1:100
%     z_all(:,:,ii)=z;
% end

disp('load data')
data_count=0;
data_objs = {};
StatsOPMs = {};
pinwheel_pos = {};
tracker_obj = pinwheel_tracker;
TrackingResults = {};
for ii = 1:N%length(data_info.exp_data)
    if ~isempty(data_info.exp_data(ii).binocular)
        %% load data
        data_count=data_count+1;
        [data_objs{data_count},~,~]=loadMouseLemurData(data_info,data_path,ii,1);
        data_objs{data_count}.prepare_samples_array(bootstrapsamples)
        disp(['loaded experiment' num2str(ii)])

        z = data_objs{data_count}.filter_map(data_objs{data_count}.read_map());
        [~,~,~,pinwheel_pos{data_count}.x,pinwheel_pos{data_count}.y,~, ~] = find_pinwheels(z,0,data_objs{data_count}.ROI,0);
        disp(['finished pinwheel finding experiment' num2str(ii)])

        %% track pinwheels between maps
        if data_count>2
            z1 = data_objs{data_count-1}.filter_map(data_objs{data_count-1}.read_map());
            z2 = data_objs{data_count}.filter_map(data_objs{data_count}.read_map());
            TrackingResults{data_count-1} = tracker_obj.interpolate(z1,z2,data_objs{data_count-1}.ROI.*data_objs{data_count}.ROI);
        
%             if data_count ==2
%                 TrackingList = zeros(size(TrackingResults{data_count-1}.ini,1),length(data_objs));
%                 TrackingList(:,1:2)=TrackingResults{data_count-1}.ini;
%             else
%                 IDs = TrackingList(:,data_count-1);
%                 for jj = 1:length(IDs)
%                     if IDs(jj)~=0
%                         TrackingList(jj,data_count) = TrackingResults{data_count-1}.ini(IDs(jj),2);
%                     end
%                 end 
%             end
%             disp(['tracked pinwheels between experiments' num2str(ii-1) 'and' num2str(ii)])
%             data_objs{data_count-1}=[];
        end
        
        

%         %% calc OPM and Pw Stats 
% 
%         StatsOPMs{data_count} = getStatsOPM(data_objs{data_count});
%         disp(['calculated stats for experiment' num2str(ii)])

    else 
        disp(['No binocular data in experiment ' num2str(ii)])
    end
end
disp('-----------------------')

DataFile = [ResultDataFolder 'AllDataObjsBiocular' animal num2str(SpecimenNumber) '.mat'];
save(DataFile,'data_objs','data_info','pinwheel_pos','StatsOPMs','TrackingList','TrackingResults')
load (DataFile,'data_objs','data_info','pinwheel_pos','StatsOPMs','TrackingList','TrackingResults')


%% track between mean maps
% tracker_obj = pinwheel_tracker;
% TrackingResults = {};
% for ii = 2:length(data_objs)
%     z1 = data_objs{ii-1}.filter_map(data_objs{ii-1}.read_map());
%     z2 = data_objs{ii}.filter_map(data_objs{ii}.read_map());
%     TrackingResults{ii-1} = tracker_obj.interpolate(z1,z2,data_objs{ii-1}.ROI.*data_objs{ii}.ROI);
% 
%     if ii ==2
%         TrackingList = zeros(size(TrackingResults{ii-1}.ini,1),length(data_objs));
%         TrackingList(:,1:2)=TrackingResults{ii-1}.ini;
%     else
%         IDs = TrackingList(:,ii-1);
%         for jj = 1:length(IDs)
%             if IDs(jj)~=0
%                 TrackingList(jj,ii) = TrackingResults{ii-1}.ini(IDs(jj),2);
%             end
%         end 
%     end
% end

%% calculate pinwheel distances
StepDistances = zeros(size(TrackingList),'like',0.5i);
FullDistances = zeros(size(TrackingList),'like',0.5i);
for ii = 1:size(TrackingList,2)-1
    IDs = TrackingList(:,ii);
    for jj = 1:length(IDs)
        StepDistances(jj,ii+1) = getPinwheelDistance(ii,ii+1,jj,TrackingList,pinwheel_pos);
        FullDistances(jj,ii+1) = getPinwheelDistance(1,ii+1,jj,TrackingList,pinwheel_pos);
    end
end

%% get pinwheel stats

PwProb = zeros(size(TrackingList,2));
PwCA = zeros(size(TrackingList,2));
for ii = 1:size(TrackingList,2)
    IDs = TrackingList(:,ii);
    SizesCI = getConfidenceRegionPw(StatsOPMs{ii}.pinwheel_stats,data_info.field_size_pix,0.95,false);
    for jj = 1:length(IDs)
        if IDs(jj) ~=0
            PwProb(jj,ii) = StatsOPMs{ii}.pinwheel_stats.probability(IDs(jj));
            PwCA(jj,ii) = SizesCI(IDs(jj));
        else
            PwProb(jj,ii) = nan;
            PwCA(jj,ii) = nan;
        end
    end
end

%% make uniform tracking list

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


function PwDistance = getPinwheelDistance(i1,i2,jj,TrackingList,pinwheel_pos)
    
    if (TrackingList(jj,i2)*TrackingList(jj,i1))==0
        PwDistance = nan;
    else
    x1 = pinwheel_pos{i1}.x(TrackingList(jj,i1));
    x2 = pinwheel_pos{i2}.x(TrackingList(jj,i2));
    y1 = pinwheel_pos{i1}.y(TrackingList(jj,i1));
    y2 = pinwheel_pos{i2}.y(TrackingList(jj,i2));
    

    %PwDistance = sqrt((x2-x1)^2+(y2-y1));
    PwDistance = (x2-x1)+1i*(y2-y1);
    end
end

function [data_obj,data,BloodVesselImg]=loadMouseLemurData(data_info,data_path,day,experiment)
    
    %% input
    animal = 'mouse lemur';
    if nargin <3
        day = 1;
    end
    if nargin<4
        experiment = 1;
    end
    
%     %% get data info
%     [data_info,data_path] = info_handle(animal,SpecimenNumber,AnimalDataFolder);
%     if isfield(data_info,'pix_per_mm')
%         data_info.pixels_per_mm = data_info.pix_per_mm;
%     else
%         data_info.pix_per_mm = data_info.pixels_per_mm;
%     end

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

