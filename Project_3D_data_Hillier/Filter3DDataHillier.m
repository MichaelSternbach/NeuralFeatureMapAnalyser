clear all
close all
disp('Start OPM analysis')
disp('-------------------------------------------')
%% add path to OPM processing functions 
addSuperdirSubfolder('OPM_Processing')


%% animal parameter
animal = 'Cat_PwTracking';


%% general parameter data evaluation
getCI = false;
bootstrapsamples = 100;


%% data folder
AnimalDataFolder = '/home/michael/Cloud/Cloud/PhD/data/CatHillier/3D_Data/';
ResultDataFolder = ['Data/' lower(animal) '/'];
mkdir(ResultDataFolder)

%% plot parameter
depth = 1;
FigureFileName = [ResultDataFolder 'ResultOPM_' animal '.eps'];
border_ColumnSpacing = 20;
position_ColumnSpacing_txt = [10 10];
x_hypercolumns = 1;


%% load data
disp('load data1')

FileName = 'Leo_20250408_phase_maps_params_masks.mat';
%FileName='Lulu_20241030_phase_maps_params_masks.mat';%Lulu1
%FileName = 'Lulu_20241211_phase_maps_params_masks.mat';%Lulu2
% 'Lulu_20240711_phase_maps_params_masks.mat';
pixel_size_um = [200 60 125];
pixel_size_new_um = min(pixel_size_um);
FileData = load([AnimalDataFolder FileName]);
disp('-----------------------')
FigureFolder = 'FiguresTracking3D/Lulu/';

%%Lulu1
%load('3DTrackingData.mat','PinwheelTracks_p','PinwheelTracks',"ROI","ROI_p")
%%Lulu2
% load('3DTrackingData2.mat','PinwheelTracks_p','PinwheelTracks',"ROI","ROI_p")

%% make data
data3D = makeOPMFromPhaseMap3D(FileData.per_cycle_phase_maps);

% data3D_uniform = resample_to_uniform_density(data3D, pixel_size);

%% make map and ROI
z =mean(data3D,4);
ROI = FileData.volume_mask==1;

%% plot orig data
ROI2D = reshape(ROI(depth,:,:),size(ROI,[2 3]));
figure;
plot_map(reshape(z(depth,:,:),size(z,[2 3])),ROI2D,0,1)
title('unfiltered')

% ROI2D = reshape(prod(ROI,1),size(ROI,[2 3]));
% figure; plot_map(ROI2D)

%% make pixel density uniform
[z, ROI] = resample_to_uniform_density(z, pixel_size_um,ROI,pixel_size_new_um);
ROI= (ROI==1);


%% plot uniform data
ROI2D = reshape(ROI(depth,:,:),size(ROI,[2 3]));
figure;
plot_map(reshape(z(depth,:,:),size(z,[2 3])),ROI2D,0,1)
title('unfiltered uniform')


%% prepare data_info
data_info.animal = animal;
data_info.ID = 'Test1';
data_info.stim_order = FileData.orientations_in_order;%[0,22.500000000000000,45,67.500000000000000,90,1.125000000000000e+02,135,1.575000000000000e+02];
data_info.pix_per_mm = 1000/pixel_size_new_um;%16.666666666666670;
data_info.pixels_per_mm = data_info.pix_per_mm;


%% calculate powerprofile
% profile_range_mm = [0.1 2];
% profile_step_mm = 0.02;
% 
% power_profile = define_filter_settings3D(data_info,ROI,z,profile_range_mm,profile_step_mm);
% 
% figure
% plot(power_profile.scale_mm,power_profile.values);
% xlabel('Scale in mm')
% ylabel('Power')
% xlim([0 max(profile_range_mm)])
% set(gca,'fontsize',15)

lowpass_cutoff = 0.35;
highpass_cutoff = 1.25;


%% filter 3D data
rise = 0.05;
z_filtered = filter_map3D(z,ROI, data_info.pixels_per_mm , highpass_cutoff, lowpass_cutoff, rise);


z3D = reshape(z_filtered(depth,:,:),size(z_filtered,[2 3]));
z3D=(z3D-mean(z3D(ROI2D)))/std(z3D(ROI2D));
figure;
plot_map(z3D,ROI2D)
title('filtered in 3D')



%% Track pinwheels
% tic
% % layers = 10:130;
layers = 10:size(z_filtered,1);
tracker_obj= pinwheel_tracker;
PinwheelTracks = trackPws3D(z_filtered(layers,:,:),ROI(layers,:,:),tracker_obj);
% toc

% save('3DTrackingData.mat','PinwheelTracks',"ROI")

%% Plot layers
% Xlim=[123 152];
% Ylim = [38 62];
% Xlim=[171 200];
% Ylim = [4 29];
% for depth = 1:size(ROI,1)
%     plotLayer(z_filtered,ROI,depth,PinwheelTracks)
%     xlim(Xlim)
%     ylim(Ylim)
%     savefig([FigureFolder 'depth' num2str(depth) 'OPM.fig' ])
%     close all
% end

%% plot Tracks and calc stats
figure;
TrackStats = plotPwTracks3D(PinwheelTracks,5,ROI,false,false);

figure;
TrackStats = plotPwTracks3D(PinwheelTracks,5,ROI,false,true);
% 
% xlim(Xlim)
% ylim(Ylim)
% 
% % Remove text outside axis limits
% ax = gca; % Get current axis handle
% xLimits = xlim(ax);
% yLimits = ylim(ax);
% zLimits = zlim(ax);
% 
% % Find and remove text objects outside the limits
% allText = findall(ax, 'Type', 'text'); % Find all text objects
% for i = 1:length(allText)
%     pos = allText(i).Position; % Get text position
%     if pos(1) < xLimits(1) || pos(1) > xLimits(2) || ...
%        pos(2) < yLimits(1) || pos(2) > yLimits(2) || ...
%        pos(3) < zLimits(1) || pos(3) > zLimits(2)
%         delete(allText(i)); % Remove text outside limits
%     end
% end
% 
% hold off;
% zlabel('depth')
% savefig([FigureFolder 'PwTracks3D.fig' ])

%view([1, 0, 0])

%% Plot Hitsograms
figure;
histogram(TrackStats.TrackSteps(TrackStats.TrackSteps>5))
xlabel('TrackSteps')

figure;
histogram(TrackStats.Start)
xlabel('Start')

figure;
histogram(TrackStats.End)
xlabel('End')

figure;
histogram(TrackStats.TrackLength./data_info.pix_per_mm)
xlabel('TrackLength [mm]')


figure;
histogram(mod(TrackStats.angle2/(2*pi)*360,360))
xlabel('angle [°]')
xlim([0 100])




%% permute data
change = [2 1 3];
z_filtered_p = permute(z_filtered,change);
ROI_p = permute(ROI,change);

% plotLayer(z_filtered_p,ROI_p,120)

%% Track pinwheels
% tracker_obj= pinwheel_tracker;
% tic
% %layers = 10:130;%Lulu1
% layers = 1:120;%Lulu2
% PinwheelTracks_p = trackPws3D(z_filtered_p(layers,:,:),ROI_p(layers,:,:),tracker_obj);
% toc

% save('3DTrackingData.mat','PinwheelTracks_p','PinwheelTracks',"ROI","ROI_p")

figure;
TrackStats_p = plotPwTracks3D(PinwheelTracks_p,5,ROI_p,false,false);
figure;
TrackStats_p = plotPwTracks3D(PinwheelTracks_p,5,ROI_p,false,true);

% Tracks_p = switchZandX(TrackStats_p.Tracks);


%% Plot Hitsograms
figure;
histogram(TrackStats_p.TrackSteps(TrackStats_p.TrackSteps>5))
xlabel('TrackSteps')

figure;
histogram(TrackStats_p.Start)
xlabel('Start')

figure;
histogram(TrackStats_p.End)
xlabel('End')

figure;
histogram(TrackStats_p.TrackLength./data_info.pix_per_mm)
xlabel('TrackLength [mm]')


figure;
histogram(mod(TrackStats_p.angle2/(2*pi)*360,360))
xlabel('angle [°]')
xlim([0 100])

%% compare trajectories

[pairing, mean_distance, directions_set1, directions_set2] = analyzeTrajectories(TrackStats.Tracks, Tracks_p, 5);

% 
% load('3DTrackingData.mat','PinwheelTracks_p','PinwheelTracks',"ROI","ROI_p")
% 
% %% Compare different tracking
% close all
% a=0;
% figure;
% %TrackComparison = ComparePwTracks3D(PinwheelTracks,PinwheelTracks_p,5,ROI,min(layers));
% %PinwheelTracks,PinwheelTracksPermutated,minSteps,ROI,base_permutated
% TrackComparison = ComparePwTracks3D(PinwheelTracks,5,ROI,a);
% figure;
% TrackStats = plotPwTracks3D(PinwheelTracks,5,ROI);

% vercoplotLayers(z_filtered(layers,:,:),ROI(layers,:,:),PinwheelTracks)
% 
% 
% %% plot 3D contours
% AngleMap = angle(z_filtered);
% plot3DContour(AngleMap, pi,ROI);
% 
% %% get 2D filtered Map
% % data2D = makeDataFromPhaseMap(reshape(FileData1.per_cycle_phase_maps(:,depth,:,:), ...
% %     size(FileData1.per_cycle_phase_maps,[1 3 4])),data_info.stim_order);
% data2D = ones([size(z_filtered,[2 3]) length(data_info.stim_order) 2]);
% data_obj = data_handle_corrected(data_info,data2D);
% data_obj.set_filter_parameters('lowpass', lowpass_cutoff)
% data_obj.set_filter_parameters('highpass', highpass_cutoff)
% z2D = data_obj.filter_map(reshape(z(depth,:,:),size(z_filtered,[2 3])));
% 
% 
% 
% 
% %% compare 3D and 2D filtering
% 
% z3D = reshape(z_filtered(depth,:,:),size(z_filtered,[2 3]));
% z3D=(z3D-mean(z3D(ROI2D)))/std(z3D(ROI2D));
% figure;
% plot_map(z3D,ROI2D)
% title('filtered in 3D')
% 
% figure;
% plot_map(z2D,ROI2D)
% title('filtered in 2D')
% 
% DifferenzeMap = abs(z2D-z3D)./abs(z2D);%/mean(abs(z2D),"all");
% 
% 
% figure;
% histogram(DifferenzeMap(ROI2D))
% xlabel('rel. difference')
% 
% DiffMax = 3;
% xlim([0 3])
% DifferenzeMap(DifferenzeMap>DiffMax)=DiffMax;
% DifferenzeMap(~ROI2D)=0;
% 
% figure;
% plot_mapAbs(DifferenzeMap)
% title('rel. difference')





%% FUnctions

function Plot2DLayer(z,ROI,PinwheelTracks,layer)
ROI2D = reshape(ROI(layer,:,:),size(ROI,[2 3]));
figure;
plot_map(reshape(z(layer,:,:),size(z,[2 3])),ROI2D,0,1)
title(['layer ' num2str(layer)])
hold on
plot(PinwheelTracks.x(:,layer),PinwheelTracks.y(:,layer),'xw')
hold on
offset = 1;

    for ii = size(PinwheelTracks.x,1)
    
        text(PinwheelTracks.x(ii,layer)+offset,PinwheelTracks.x(ii,layer)+offset,num2str(ii))
        hold on
    
    end
end


function plotLayers(z_filtered,ROI,PinwheelTracks)
    for ii = 1:size(z_filtered,1)
        plotLayer(z_filtered,ROI,ii,PinwheelTracks)
    end
end

function plotLayer(z_filtered,ROI,depth,PinwheelTracks)
    ROI2D = reshape(ROI(depth,:,:),size(ROI,[2 3]));
    z3D = reshape(z_filtered(depth,:,:),size(z_filtered,[2 3]));
    z3D=(z3D-mean(z3D(ROI2D)))/std(z3D(ROI2D));
    figure;
    plot_map(z3D,ROI2D)
    
    if nargin >3
        hold on
        plot(PinwheelTracks.x(:,depth),PinwheelTracks.y(:,depth),'xw')
        
        label_offset = 1;
        for ii = 1:size(PinwheelTracks.label,1)
            label = PinwheelTracks.label(ii,depth);
            if label ~=0
                x = PinwheelTracks.x(label,depth)+label_offset;
                y = PinwheelTracks.y(label,depth)+label_offset;
    
                text(x,y,num2str(label),'Color','white')
            end
        end
    end
    title(['depth ' num2str(depth)])
end

function z_filtered = BandpassFilterOPM3D(z,lowpass_cutoff,highpass_cutoff,data_info)
    steepness = 0.05;
    
    lowpass_cutoff_pix = lowpass_cutoff * data_info.pix_per_mm;
    highpass_cutoff_pix = highpass_cutoff * data_info.pix_per_mm;

    z_lowpass = filterMap3D(z, lowpass_cutoff_pix,steepness);
    z_filtered = z_lowpass-filterMap3D(z, highpass_cutoff_pix,steepness);
end



function power_profile = define_filter_settings3D(data_info,ROI,z,profile_range_mm,profile_step_mm)
    %(experiment_num)
    % Finding appropiate filter settings is crucial in the data analysis.
    % This function calculates the radial profile of the power spectrum of the 
    % maps in each trial to use as a reference to set initial filter cut-offs. 
    % Since the cutt-offs are defined in mm, the function calculates the
    % profile scale in mm instead of k (wavenumber).
    % The results are all plotted in a figure and saved as power_profile 
    % in exp_info.mat.
    % When the map is of good quality, the power is low for small mm scales,
    % has a few large peaks corresponding to the typical scale of the map around
    % 0.6-1.0 mm, decreases afterwards and starts fluctuating depending on the
    % larger structures of the layout.
    % A good highpass is where the first peaks have dropped and before the power
    % rises again.
    % A good lowpass is where the first peak starts raising for all maps.
    % The lowpass is crucial for the pinwheel position, so this initial
    % estimate will be refined in a second step later on.
    % The defined cut-offs should be saved in make_info_files.m under settings:
    % e.g.
    % data_info.settings.lowpass_mm = 0.47;
    % data_info.settings.highpass_mm = 1.5;
    %%
    
    % parameters for profile
    % profile_range_mm = [0.1 5];
    % profile_step_mm = 0.01;
    
%     z(ROI)=0;
    z = z.*double(ROI);

    if length(profile_range_mm)>2
        profile_scale_mm = profile_range_mm;
    else
        profile_scale_mm = profile_range_mm(1):profile_step_mm:profile_range_mm(2);
    end
    
    stim_order = data_info.stim_order;
    power_profile.scale_mm = [];
    power_profile.values = [];
    
    
    %% calculate power spectrum
    
    power_profile.values = calculate_power_profile_3d(z,profile_scale_mm*data_info.pix_per_mm);
    power_profile.scale_mm = profile_scale_mm;
    
    
    
    %% normalize powerspectrum
    power = mean(abs(z).^2,'all');
    power_profile.k_mm_inv = 1./power_profile.scale_mm;
    scale_mm_fine = linspace(min(power_profile.k_mm_inv), max(power_profile.k_mm_inv), 1000);  % Create a fine grid over the range of x
    values_fine = interp1(power_profile.k_mm_inv, power_profile.values, scale_mm_fine, 'linear');  % Linear interpolation on the fine grid
    integral_power = trapz(scale_mm_fine, values_fine); 
    
    power_profile.values_kspace = power_profile.values/integral_power*power;

end
