addpath /home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/AGWolfOPMDataPipelineCode
addpath /home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/OnlineAnalysisOPM

%clear all

%% parameter

experiment_num = 1;
animal = 'dunnart';

min_lowpass_mm = 0.3;
max_lowpass_mm = 0.8;

n_steps_lowpass= 10;


%% load data

[data_info,data_path] = info_handle(animal,experiment_num);
data_info.pixels_per_mm = data_info.pix_per_mm;

data = NoAveragePreProcessRawDataJason(data_info.expIds,data_info.refWin,data_info.sigWin,data_info.partId,data_path,data_info.ID);
data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
data_obj.set_filter_parameters('lowpass',.5)
data_obj.set_filter_parameters('highpass',1.2)


%% Make tracker 
tracker_obj = pinwheel_tracker;

%% run pinwheel tracking

[PwFilterVelocity,PwFilterVelocityMean,PWNumber,lowpass_cutoffs] = getPWFilterVelocityAndNumber(data_obj,tracker_obj,min_lowpass_mm,max_lowpass_mm,n_steps_lowpass);

figure
plot(lowpass_cutoffs(1:n_steps_lowpass-1), PwFilterVelocityMean)%(2:n_steps_lowpass)
hold on
for ii = 1:(n_steps_lowpass-1)
    y = PwFilterVelocity{ii};
    x = ones(size(y))*lowpass_cutoffs(ii);
    plot(x,y,'+')
    hold on
    
end
xlabel('lowpass cut-off [mm]')
ylabel('pinwheel filter velocity [pixels/mm]')
% 
% figure
% title('pinwheel number')
% plot(lowpass_cutoffs, PWNumber)
% xlabel('lowpass cut-off [mm]')
% ylabel('pinwheel number')

%% alternative way pw number

%filters = find_lowpassPwNumber(data_obj,data_info,data_obj.info.settings.highpass_mm,data_obj.info.settings.lowpass_mm,lowpass_cutoffs);
%plot(filters.global_plateau.lowpass_vs_density(:,1),filters.global_plateau.lowpass_vs_density(:,2))

%% check OrientationCI

[meanCI,CIs,lowpass_cutoffs] = getFilterOrientationCI(data_obj,min_lowpass_mm,max_lowpass_mm,n_steps_lowpass);

figure
plot(lowpass_cutoffs, meanCI)%(2:n_steps_lowpass)
xlabel('lowpass cut-off [mm]')
ylabel('mean Orientation Preference  CI [°]')

%% function

function [meanCI,CIs,lowpass_cutoffs] = getFilterOrientationCI(data_obj,min_lowpass_mm,max_lowpass_mm,n_steps_lowpass)
    Bootstrapsamples = 100;
    alpha=0.05;
    apply_filter = true;
    
    data_obj.prepare_samples_array(Bootstrapsamples)
    
    orig_lowpass = data_obj.filter_parameters.lowpass;
    
    %% define lowpass_cutoffs
    
    lowpass_cutoffs = linspace(min_lowpass_mm,max_lowpass_mm,n_steps_lowpass);
    
    %% prep variables
    meanCI = zeros(size(lowpass_cutoffs));
    CIs = cell(size(lowpass_cutoffs));
    
     %% make 1 maps
    data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(1))
    orientation_stats = get_orientation_stats(data_obj,alpha,apply_filter);
    meanCI(1) = getMeanCI(orientation_stats,data_obj.ROI);
    CIs{1} = orientation_stats;
    
    %% iterate through cut offs
    disp('iterate over different lowpass cut offs')
    for ii = 2:n_steps_lowpass
        disp(ii)
        data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(ii))
        orientation_stats = get_orientation_stats(data_obj,alpha,apply_filter);
        meanCI(ii) = getMeanCI(orientation_stats,data_obj.ROI);
        CIs{ii} = orientation_stats;        
    end
    
    
    %% revert to original lowpass
    data_obj.set_filter_parameters('lowpass',orig_lowpass)
end

function meanCI = getMeanCI(orientation_stats,ROI)
    CI = abs(angle(orientation_stats(:,:,3)./orientation_stats(:,:,1)))/pi*90;
    meanCI = mean(CI(ROI),1:2);

end


function [PwFilterVelocity,PwFilterVelocityMean,PWNumber,lowpass_cutoffs] = getPWFilterVelocityAndNumber(data_obj,tracker_obj,min_lowpass_mm,max_lowpass_mm,n_steps_lowpass)
    orig_lowpass = data_obj.filter_parameters.lowpass;
    
    %% define lowpass_cutoffs
    
    lowpass_cutoffs = linspace(min_lowpass_mm,max_lowpass_mm,n_steps_lowpass);
    
    
    %% prep variables
    
    PwFilterVelocityMean = zeros([1,n_steps_lowpass-1]);
    PWNumber = zeros([1,n_steps_lowpass]);
    
    %% make 1 maps
    data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(1))
    z_1 = data_obj.filter_map(data_obj.read_map());
    pinwheels_1 =tracker_obj.find_pinwheels(z_1,data_obj.ROI);
    PWNumber(1) = getPWNumer(pinwheels_1);
    
    %% iterate through cut offs
    
    for ii = 2:n_steps_lowpass
        
        %% make maps
        data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(ii))
        z_2 = data_obj.filter_map(data_obj.read_map());
        pinwheels_2 =tracker_obj.find_pinwheels(z_2,data_obj.ROI);
        PWNumber(ii) = getPWNumer(pinwheels_2);
        
        %% track pinwheels
        tracking = tracker_obj.interpolate(z_1,z_2,data_obj.ROI);
        
        %% calclate PW velocity
        PwFilterVelocity{ii-1} = getPwPositionDifference(tracking,pinwheels_1,pinwheels_2)/abs(lowpass_cutoffs(ii)-lowpass_cutoffs(ii-1));
        PwFilterVelocityMean(ii-1) = mean(PwFilterVelocity{ii-1},'all');
        
        z_1 = z_2;
        pinwheels_1 =pinwheels_2;
        

    end
    
    %% revert to original lowpass
    data_obj.set_filter_parameters('lowpass',orig_lowpass)
end

function PWDistances = getPwPositionDifference(tracking,pinwheels_1,pinwheels_2)
    PWNumber = getPWNumer(pinwheels_1);
    PWDistances = [];
    for ii = 1:PWNumber
       pw1 = tracking.ini(ii,1);
       x1 =  pinwheels_1.x(pw1);
       y1 =  pinwheels_1.y(pw1);
       
       pw2 = tracking.ini(pw1,2);
       if pw2 ~= 0
           x2 =  pinwheels_2.x(pw2);
           y2 =  pinwheels_2.y(pw2);

           PWDistances = [PWDistances,sqrt((x1-x2)^2+(y1-y2)^2)];
       end
    end
    %MeanPwPositionDifference = mean(PWDistances,'all');
end

function PWNumber = getPWNumer(pinwheels)
    PWNumber = size(pinwheels.label,1);
end