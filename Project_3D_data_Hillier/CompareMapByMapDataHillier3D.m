clear all
close all
disp('Start OPM analysis')
disp('-------------------------------------------')
%% add path to OPM processing functions 
addpath /home/michael/Cloud/Cloud/PhD/MarsupialData/OrientationPrefernceMapProcessing/OPM_Processing

%% data folder
AnimalDataFolder = '/home/michael/Cloud/Cloud/PhD/data/CatHillier/3D_Data/';

%% load data
disp('load data1')

FileNamePhaseMap='Lulu_20250218_phase_maps_params_masks.mat';
FileNameAzimuthMap='Lulu_20250218_azim_retinotopy_deg_maps_corrected.mat';
FileNameElevationMap = 'Lulu_20250218_elev_retinotopy_deg_maps_corrected.mat';
pixel_size_um = [200 60 125];
pixel_size_new_um = min(pixel_size_um);
OPMData = load([AnimalDataFolder FileNamePhaseMap]);
AzimutData = load([AnimalDataFolder FileNameAzimuthMap]);
ElevationData = load([AnimalDataFolder FileNameElevationMap]);

AzimutData = AzimutData.retinotopy_map;
ElevationData = ElevationData.retinotopy_map;

disp('-----------------------')


%% make data
data3D = makeOPMFromPhaseMap3D(OPMData.per_cycle_phase_maps);

%% make map and ROI
z3D =mean(data3D,4);
ROI3D = OPMData.volume_mask==1;



%% prepare data_info
data_info.animal = 'cat';
data_info.ID = 'Test1';
data_info.stim_order = OPMData.orientations_in_order;%[0,22.500000000000000,45,67.500000000000000,90,1.125000000000000e+02,135,1.575000000000000e+02];
data_info.pix_per_mm = 1000/pixel_size_new_um;%16.666666666666670;
data_info.pixels_per_mm = data_info.pix_per_mm;

data_info.settings.lowpass_mm = 0.35;
data_info.settings.highpass_mm = 1.25;


%% plot 1 layer
evaluateLayer(1,z3D,ROI3D,AzimutData,ElevationData,data_info,pixel_size_um,pixel_size_new_um)

evaluateLayer(1,z3D,ROI3D,AzimutData,ElevationData,data_info,pixel_size_um,pixel_size_new_um,true)

evaluateLayer(5,z3D,ROI3D,AzimutData,ElevationData,data_info,pixel_size_um,pixel_size_new_um)


evaluateLayer(10,z3D,ROI3D,AzimutData,ElevationData,data_info,pixel_size_um,pixel_size_new_um)



%% FUnctions

function evaluateLayer(layer,z3D,ROI3D,AzimutData,ElevationData,data_info,pixel_size_um,pixel_size_new_um,filterVF)

    if nargin <9
        filterVF = false;
    end

    z = z3D(layer,:,:);
    z = reshape(z,size(z3D,[2 3]));
    
    ROI = ROI3D(layer,:,:);
    ROI = reshape(ROI,size(ROI3D,[2 3]));
    
    AzimutData = AzimutData(layer,:,:);
    AzimutData = reshape(AzimutData,size(AzimutData,[2 3]));
    
    ElevationData = ElevationData(layer,:,:);
    ElevationData = reshape(ElevationData,size(ElevationData,[2 3]));
    
    %% make pixel density uniform
    ROI_old = ROI;
    [z, ROI] = resample_to_uniform_density(z, pixel_size_um(2:3),ROI,pixel_size_new_um);
    ROI= (ROI==1);
    
    [AzimutData, ~] = resample_to_uniform_density(AzimutData, pixel_size_um(2:3),ROI_old,pixel_size_new_um);
    [ElevationData, ~] = resample_to_uniform_density(ElevationData, pixel_size_um(2:3),ROI_old,pixel_size_new_um);
    
    
    data2D = ones([size(z,[1 2]) length(data_info.stim_order) 2]);
    data_obj = data_handle_corrected(data_info,data2D,ROI);
    % data_obj.set_filter_parameters('lowpass', lowpass_cutoff)
    % data_obj.set_filter_parameters('highpass', highpass_cutoff)
    z = data_obj.filter_map(z);
    
    if filterVF
        cutoff = 2;
        width = 5;
        AzimutData = fermi_lowpass_2d(AzimutData,cutoff,width);
        ElevationData = fermi_lowpass_2d(ElevationData,cutoff,width);
    end
    
    %% get pinwheel density
%     [average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient] = get_column_spacingManuel(z,ROI,data_info.pix_per_mm,0.3,2,0.05);
    [PwInfo.NumberPw,PwInfo.aniso,PwInfo.x_angle,PwInfo.PWxList,PwInfo.PWyList,PwInfo.signList, PwInfo.contours] = find_pinwheels(z,0,ROI);
%     PwInfo.NumHypercolumns = sum(ROI,'all')/(data_obj.info.pix_per_mm*average_spacing_mm)^2;
%     
%     sigma = 0.1;
%     PwInfo.LocalPwDensityFixedFilter=getLocalPwDensityFixedFilter(data_obj,PwInfo,local_spacing_mm,sigma);
%     PwInfo.WeightedPwDensityFixedFilter = mean(PwInfo.LocalPwDensityFixedFilter(ROI));
%     
    %% prepare coordinates
            
    
    X = AzimutData;
    Y = ElevationData;
    % X_range = [min(X(ROI),[],'all') max(X(ROI),[],'all')];
    % Y_range = [min(Y(ROI),[],'all') max(Y(ROI),[],'all')];
    X_range = [-45 45];
    Y_range = [-30 30];
    
    %% get pixel density from minimum distance between points
    pixel_density = 1/(min(diff(X_range)/size(X,2),diff(Y_range)/size(Y,1)));
    
    
    z_VF = interpolate_2d_array(z(ROI), X(ROI), Y(ROI), X_range, Y_range, pixel_density);
    %ROI_VF = interpolate_2d_array(ROI2D, X, Y, X_range, Y_range, pixel_density);

%     %% get column spacing z visual field
%     [average_spacing_mm_VF,local_spacing_mm_VF,newROI,WavletCoefficient] = get_column_spacingManuel(z_VF,ROI,1,0.3,2,0.05);
%     
    %% plot results
    % figure;plot_map(z,ROI,0,1)
    % figure;plot_mapAbs(X)
    % figure;plot_mapAbs(Y)
    % figure;plot_map(z_VF)
    
    figure;
    t = tiledlayout(2,3);
    title(t,'Layer ' + string(layer))
    
    nexttile
    plot_map(z,ROI,0,1)
    hold on
    plot(PwInfo.PWxList,PwInfo.PWyList,'xw')
    title('z')
    
%     nexttile
%     % plot column spacing, pinwheel number, number hypercolumns and density as txt
%     text(0.1,0.9,['Column spacing: ' num2str(average_spacing_mm) ' mm'])
%     text(0.1,0.8,['Pinwheel number: ' num2str(PwInfo.NumberPw)])
%     text(0.1,0.7,['Number of hypercolumns: ' num2str(PwInfo.NumHypercolumns)])
%     text(0.1,0.6,['pinwheel density: ' num2str(PwInfo.NumberPw/PwInfo.NumHypercolumns) ' 1/mm^2'])
%     title('Column spacing and pinwheel density')
    
    ax = nexttile;
    %plot_mapAbs(map,Title,maxMap,minMap,ROI,ax)
    plot_mapAbs(X,'Azimut',max(X(ROI),[],'all'),min(X(ROI),[],'all'),ROI,ax)
    
    ax = nexttile;
    plot_mapAbs(Y,'Elevation',max(Y(ROI),[],'all'),min(Y(ROI),[],'all'),ROI,ax)
    
    nexttile
    plot_mapVF(z_VF)
    % add tick positions with values based on X and Y
    N_ticks = 5;
    set(gca,'XTick',linspace(1,size(z_VF,2),N_ticks))
    set(gca,'YTick',linspace(1,size(z_VF,1),N_ticks))
    set(gca,'XTickLabel',linspace(X_range(1),X_range(2),N_ticks))
    set(gca,'YTickLabel',linspace(Y_range(1),Y_range(2),N_ticks))
    set(gca,'ydir','normal')

    nexttile;
    OriScatterPlot(z(ROI),X(ROI),Y(ROI))
%     set(gca,'ydir','reverse')

%     ax = nexttile;
%     plot_mapAbs(local_spacing_mm_VF,'spacila scale [°]',max(local_spacing_mm_VF(ROI),[],'all'),min(local_spacing_mm_VF(ROI),[],'all'),ROI,ax)

    nexttile;
    % Define the bin edges
    x = X(ROI);
    y = Y(ROI)
    x_edges = linspace(min(x), max(x), 50); % 50 bins in x-direction
    y_edges = linspace(min(y), max(y), 50); % 50 bins in y-direction
    
    % Count number of points in each bin
    counts = histcounts2(x, y, x_edges, y_edges);
    
    % Plot heat map
    imagesc(x_edges, y_edges, counts.');
    axis xy; % So that y increases upwards
    colorbar;
    xlabel('X');
    ylabel('Y');

    %% plot histogram bins
    
    plotZHistograms(x, y, angle(z_VF(ROI)), 6)

end

function V_interp = interpolate_2d_array(V, X, Y, X_range, Y_range, pixel_density)
    % Determine the size of the output array
    x_min = X_range(1);
    x_max = X_range(2);
    y_min = Y_range(1);
    y_max = Y_range(2);
    
    % Generate the new grid based on the given pixel density
    xq = linspace(x_min, x_max, round((x_max - x_min) * pixel_density));
    yq = linspace(y_min, y_max, round((y_max - y_min) * pixel_density));
    [Xq, Yq] = meshgrid(xq, yq);
    
    % Perform linear interpolation
    V_interp = griddata(X, Y, V, Xq, Yq, 'linear');
end

function plotZHistograms(x, y, z, n_bins)
% plotZHistograms Bins x and y data and plots histograms of z for each bin.
%
% Usage:
%   plotZHistograms(x, y, z)
%   plotZHistograms(x, y, z, n_bins)
%
% Inputs:
%   x, y, z  - Vectors of the same length.
%   n_bins   - (Optional) Scalar specifying the number of bins in x and y.
%              Default is 5.

    if nargin < 4
        n_bins = 5; % Default number of bins if not provided
    end

    % Define bin edges
    xEdges = linspace(min(x), max(x), n_bins+1);
    yEdges = linspace(min(y), max(y), n_bins+1);

    % Bin indices for each data point
    [~,~,xBin] = histcounts(x, xEdges);
    [~,~,yBin] = histcounts(y, yEdges);

    % Set up figure
    figure;
    plotIdx = 1;

    for i = 1:n_bins
        for j = 1:n_bins
            % Find indices of points in the current (x,y) bin
            inBin = (xBin == i) & (yBin == j);

            % Get corresponding z values
            z_in_bin = z(inBin);

            % Get the x and y range for the current bin
            xRange = sprintf('[%.2f, %.2f]', xEdges(i), xEdges(i+1));
            yRange = sprintf('[%.2f, %.2f]', yEdges(j), yEdges(j+1));

            % Create subplot
            subplot(n_bins, n_bins, plotIdx);
            
            if ~isempty(z_in_bin)
                histogram(z_in_bin, 'Normalization', 'count');  % Absolute counts
                xlabel('z');
                ylabel('Count');
            else
                axis off; % Turn off axis if no data
            end
            
            % Title with x and y ranges
            title(sprintf('x: %s, y: %s', xRange, yRange));
            plotIdx = plotIdx + 1;
        end
    end

    sgtitle('Histograms of z per (x,y) bin');

end


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

function filtered_img = fermi_lowpass_2d(img, cutoff, width)
    % FERMI_LOWPASS_2D Applies a 2D Fermi low-pass filter to an image.
    % 
    % filtered_img = fermi_lowpass_2d(img, cutoff, width)
    % img    - Input 2D image or array
    % cutoff - Cutoff frequency (in pixels, relative to image size)
    % width  - Transition width of the filter (in pixels)
    %
    % Returns:
    % filtered_img - The filtered image in the spatial domain

    % Get image size
    [rows, cols] = size(img);
    
    % Create frequency grid
    [u, v] = meshgrid(-floor(cols/2):floor((cols-1)/2), -floor(rows/2):floor((rows-1)/2));
    freq_radius = sqrt(u.^2 + v.^2);
    
    % Define Fermi filter in frequency domain
    fermi_filter = 1 ./ (1 + exp((freq_radius - cutoff) / width));
    
    % Apply FFT to image and shift to centre
    img_fft = fftshift(fft2(img));
    
    % Apply filter in frequency domain
    filtered_fft = img_fft .* fermi_filter;
    
    % Inverse FFT to get the filtered image
    filtered_img = real(ifft2(ifftshift(filtered_fft)));
end

function OriScatterPlot(Z,X,Y)
    cm = makeColormap('orientation',16);   
    % normalize from 0-1
    ori=(angle(Z))/(2*pi);ori(ori<0)=1+ori(ori<0);
    ori=ori(:);
    
    % find to which color interval it belongs
    intervals=linspace(0,1,2*size(cm,1)+1);
    [~,ori] = histc(ori,intervals(2:2:end));
    ori = ori+1;
    ori(ori==size(cm,1)+1)=1;

    % Plot with RGB colours
    scatter(X, Y, 36, ori, 'filled');
end

function varargout=plot_mapVF(data,ROI,ref_sel,black_roi,smoothing,color_type)
    % [h,fig2plot]=plot_map(data,ROI,ref_sel,black_roi,smoothing,color_type)
    % Function to make a color and intensity colored orientation preference
    % display
    % input = polar data, Region of interest
    % mapPolar(data2D,ROI)
    % created by chepe@nld.ds.mpg.de
    
    % Other functions used:
    % - makeColormap
    % - smooth_pf, negative_smoothing, imresize
    
    
    %% ---------------------  read parameters
    
    data(isnan(data)) = 0;
    
    % if empty ROI, use all image
    if ~exist('ROI','var') || isempty(ROI)
        ROI=true(size(data));
    end
    
    % determine polar scaling for selectivity
    if ~exist('ref_sel','var')
        ref_sel=3*sqrt(mean(abs(data(ROI)).^2));
    elseif isempty(ref_sel)
        ref_sel=3*sqrt(mean(abs(data(ROI)).^2));
    end
    
    if ~exist('black_roi','var')
        black_roi = 1;
    end
    
    % increase size of display for exporting figures
    if exist('smoothing','var') && ~isempty(smoothing)
        if smoothing>0
            for times=1:smoothing
                data=smooth_pf(data);
            end
        elseif smoothing<0
            for times=1:abs(smoothing)
                data = negative_smoothing(data);
            end
        end
        ROI = imresize(ROI,2^abs(smoothing));
    end
    
    %-  make color maps
    if ~exist('color_type','var')
        color_type = 'std';
    end
    switch color_type
        case 'std'
            cm = makeColormap('orientation',16);        
        case 'interp'
            cm = makeColormap('orientation',24);
        case 'circ'
            cm = makeColormap('circular',50);        
    end
    mm = makeColormap('selectivity',64);
    %% --------------------------- Asign color values
    
    % -- ORIENTATION
    
    % normalize from 0-1
    ori=(angle(data))/(2*pi);ori(ori<0)=1+ori(ori<0);
    ori=ori(:);
    
    % find to which color interval it belongs
    intervals=linspace(0,1,2*size(cm,1)+1);
    [~,ori] = histc(ori,intervals(2:2:end));
    ori = ori+1;
    ori(ori==size(cm,1)+1)=1;
    
    % assign
    anglePlot = (cm(ori(:),:));
    
    % -- SELECTIVITY
    
    if ref_sel==0
        
        % if polar plot is not required
        magnitudePlot = 0.9*ones(size(anglePlot));
        
    else
        
        % threshold
        sel=abs(data)/ref_sel;
        sel(sel>1)=1;
        
        % find intensity match
        sel=ceil(size(mm,1)*sel(:));
        sel(sel==0)=1;
        
        % assign
        magnitudePlot = (mm(sel(:),:));
    end
    
    %% --------------------------- Combine
    
    if ismatrix(data)
        magnitudePlot = reshape(magnitudePlot, [size(data,1),size(data,2),3]);
        anglePlot = reshape(anglePlot, [size(data,1),size(data,2),3]);
        fig2plot=anglePlot.*magnitudePlot;
        %fig2plot=repmat(double(ROI),[1 1 3]).*anglePlot.*magnitudePlot;
        if black_roi == 1
            % black ROI
            fig2plot(repmat(~ROI,[1 1 3])) = 0;
        elseif black_roi == -1
            % white ROI
            fig2plot(repmat(~ROI,[1 1 3])) = 1;
        end
    else
        % in case of volume image, return vector
        fig2plot = magnitudePlot.*anglePlot;
    end
    
    % if only the colormap is required, return without plotting
    if nargout==2 || ~ismatrix(data)
        varargout{1} = [];
        varargout{2} = fig2plot;   
        return
    end
    
    %% --------------------------- Make figure
    
    h = image(fig2plot);
    %set(gca,'ydir','reverse')
    axis image
    
    % set(gca,'xtick',[])
    % set(gca,'ytick',[])
    
    % set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
    
    if nargout==1
        varargout{1} = h;    
    end
    
    end
    