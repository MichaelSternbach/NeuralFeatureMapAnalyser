% close all
AnimalDataFolder = '~/Cloud/';
GIF_SN_TH = 4;
FigureFolder = '~/Cloud/Cloud/PhD/MarsupialData/marsupial-data/CompareMapScale/';
size_mm = 1;

if ~ isfolder(FigureFolder)
    mkdir(FigureFolder)
end

[figCat,data_infoCat] = getMapFigure('wallaby',6,AnimalDataFolder,GIF_SN_TH, [size_mm size_mm],FigureFolder);

[figCat,data_infoCat] = getMapFigure('cat',2,AnimalDataFolder,GIF_SN_TH, [size_mm size_mm],FigureFolder);

[figCat,data_infoCat] = getMapFigure('ferret',2,AnimalDataFolder,GIF_SN_TH, [size_mm size_mm],FigureFolder);

[figCat,data_infoCat] = getMapFigure('dunnart',4,AnimalDataFolder,GIF_SN_TH, [size_mm size_mm],FigureFolder);








function [fig,data_info] = getMapFigure(animal,experiment_num,AnimalDataFolder,GIF_SN_TH, dimensions_mm,FigureFolder)
    %% animal data
    [data_info,data_path,data_obj,data,~] = getAnimalData(animal,experiment_num,AnimalDataFolder);
    if GIF_SN_TH>0
        data_obj.activateGIF(true,GIF_SN_TH)
    end
    
    
    
    %% make map borders ROI
    ROI =data_obj.ROI;
    [YROI,XROI] = find(ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);
    
    %% get center ROI
    centre = [mean([Xmin Xmax]) mean([Ymin Ymax])];

    %% calc limits
    limits = compute_plot_limits(centre, dimensions_mm, data_info.pix_per_mm);

    %% get OPM map
    z = data_obj.filter_map(data_obj.read_map());

    %% plot figure;
    fig = figure();
    plot_map(z,data_obj.ROI)
    hold on
    d = 5;
    x0 = limits.x(1) + d;
    y0 = limits.y(1) + d;
    plot([x0 x0+data_info.pix_per_mm],[y0 y0],'-w')
    title(animal)
    xlim(limits.x)
    ylim(limits.y)

    %% save as eps
    FileName = [FigureFolder '/' animal num2str(experiment_num) '_Map.eps'];
    saveas(fig,FileName,'epsc')


end

function limits = compute_plot_limits(centre, dimensions_mm, pixels_per_mm)
    % Compute pixel dimensions
    width_px = dimensions_mm(1) * pixels_per_mm;
    height_px = dimensions_mm(2) * pixels_per_mm;
    
    % Compute limits in mm
    x_min = centre(1) - width_px/ 2;
    x_max = centre(1) + width_px/2;
    y_min = centre(2) - height_px/2;
    y_max = centre(2) + height_px/2;
    
    % Return limits as a struct
    limits = struct('x', [x_min, x_max], 'y', [y_min, y_max]);
end



function combine_scaled_figures(fig_handles, pixels_per_mm_list, centre_points, display_size_mm, scale)
    % Validate input sizes
    num_figs = length(fig_handles);
    assert(length(pixels_per_mm_list) == num_figs, 'Mismatch in number of figures and pixels_per_mm values');
    assert(size(centre_points, 1) == num_figs && size(centre_points, 2) == 2, 'Invalid centre points array');
    assert(length(display_size_mm) == 2, 'Display size must have two elements');
    
    % Create a new figure with tiled layout
    tiled_fig = figure;
    tiledlayout(ceil(sqrt(num_figs)), ceil(sqrt(num_figs))); % Arrange figures in a roughly square grid
    
    % Iterate through figures and process each one
    for i = 1:num_figs
        % Extract relevant parameters
        fig = fig_handles(i);
        pixels_per_mm = pixels_per_mm_list(i);
        centre = centre_points(i, :);
        
        % Compute pixel dimensions for cropping
        crop_width_px = round(display_size_mm(1) * pixels_per_mm);
        crop_height_px = round(display_size_mm(2) * pixels_per_mm);
        
        % Compute cropping region in pixels
        x_min = round(centre(1) * pixels_per_mm - crop_width_px / 2);
        y_min = round(centre(2) * pixels_per_mm - crop_height_px / 2);
        
        % Copy and resize axes from the figure
        ax = findall(fig, 'type', 'axes');
        new_ax = nexttile;
        copyobj(ax, new_ax);
        
        % Set the limits to crop the figure properly
        xlim(new_ax, [x_min, x_min + crop_width_px] / pixels_per_mm);
        ylim(new_ax, [y_min, y_min + crop_height_px] / pixels_per_mm);
        
        % Adjust appearance
        axis(new_ax, 'equal');
    end
    
    % Adjust overall figure size based on scale factor
    if scale
        set(tiled_fig, 'Position', [100, 100, display_size_mm(1) * num_figs, display_size_mm(2) * num_figs]);
    end
end


