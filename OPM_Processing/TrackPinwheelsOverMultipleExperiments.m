function TrackPinwheelsOverMultipleExperiments(SpecimenNumber,AnimalDataFolder,ResultDataFolder,animal,N_experiments,bootstrapsamples)
    
    %% prepare parameter
    if ischar(SpecimenNumber)
        SpecimenNumber = str2num(SpecimenNumber);
    end
    % SpecimenNumber = 3;
    if nargin <4
        animal = 'mouse lemur';
    end
    
    if nargin < 3
        ResultDataFolder = ['ResultData/' lower(animal) '/' lower(animal) num2str(SpecimenNumber) '/'];
    end
    
    if nargin<6
        bootstrapsamples = 100;
    else
        if ischar(bootstrapsamples)
            bootstrapsamples = str2num(bootstrapsamples);
        end
    end
    
    
    % AnimalDataFolder = '~/CIDBN/';
    
    mkdir(ResultDataFolder)
    
    
    
    %% get data info
    [data_info,data_path] = info_handle(animal,SpecimenNumber,AnimalDataFolder);
    if isfield(data_info,'pix_per_mm')
        data_info.pixels_per_mm = data_info.pix_per_mm;
    else
        data_info.pix_per_mm = data_info.pixels_per_mm;
    end
    
    if nargin<5
        N_experiments=length(data_info.exp_data);
    else
        if ischar(N_experiments)
            N_experiments = str2num(bootstrapsamples);
        end
    end
    
    
    
    disp('load data')
    data_count=0;
    data_objs = {};
    StatsOPMs = {};
    pinwheel_pos = {};
    tracker_obj = pinwheel_tracker;
    TrackingResults = {};
    
    for ii = 1:N_experiments%length(data_info.exp_data)
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
            if data_count>1
                z1 = data_objs{data_count-1}.filter_map(data_objs{data_count-1}.read_map());
                z2 = data_objs{data_count}.filter_map(data_objs{data_count}.read_map());
                TrackingResults{data_count-1} = tracker_obj.interpolate(z1,z2,data_objs{data_count-1}.ROI.*data_objs{data_count}.ROI);
            
                if data_count ==2
                    TrackingList = zeros(size(TrackingResults{data_count-1}.ini,1),length(data_objs));
                    TrackingList(:,1:2)=TrackingResults{data_count-1}.ini;
                elseif data_count > 2 
                    IDs = TrackingList(:,data_count-1);
                    for jj = 1:length(IDs)
                        if IDs(jj)~=0
                            TrackingList(jj,data_count) = TrackingResults{data_count-1}.ini(IDs(jj),2);
                        end
                    end 
                end
                disp(['tracked pinwheels between experiments' num2str(ii-1) 'and' num2str(ii)])
                %data_objs{data_count-1}=[];
            end
            
            
    
            %% calc OPM and Pw Stats 
    
            StatsOPMs{data_count} = getStatsOPM(data_objs{data_count});
            disp(['calculated stats for experiment' num2str(ii)])
    
        else 
            disp(['No binocular data in experiment ' num2str(ii)])
        end
    end
    disp('-----------------------')
    % 
%     DataFile = [ResultDataFolder 'AllDataObjsBiocular' animal num2str(SpecimenNumber) '.mat'];
%     save(DataFile,'data_objs','data_info','pinwheel_pos','StatsOPMs','TrackingList','TrackingResults')
%     load (DataFile,'data_objs','data_info','pinwheel_pos','StatsOPMs','TrackingList','TrackingResults')
    
    
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


    %% Save result data

    DataFile = [ResultDataFolder 'TrackingResultData' animal num2str(SpecimenNumber) '.mat'];
    save(DataFile,'data_info','pinwheel_pos','StatsOPMs','TrackingList','TrackingResults')


    DataFile = [ResultDataFolder 'AllDataObjsBiocular' animal num2str(SpecimenNumber) '.mat'];
    save(DataFile,'data_objs','data_info')
    

end


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