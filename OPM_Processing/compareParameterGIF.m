function compareParameterGIF(animal,specimen_num,AnimalDataFolder,FigureFolder,ElectrodePosition,SN_th_range)
    %compareParameterGIF('dunnart',7,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/',[92 110])
    %close all; compareParameterGIF('dunnart',6,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/',[67 96])
    if nargin < 5
        ElectrodePosition = [];
    end
    if nargin <6
        SN_th_range = 0.5:0.5:4;
    end

    SizeSquare = 5;

    %% Load data
    disp('Load data...')
    [data_info, data_path, data_obj, data, BloodVesselImg] = getAnimalData(animal, specimen_num, AnimalDataFolder);
    disp('Data loaded successfully.')
    disp('-----------------------')

    %% calc plot limits
    [x_range,y_range] = getRangeXY_ROI(data_obj.ROI);

    %% set Figure file name
    FigureFileName = [FigureFolder 'ComparisonParameterGIF_' animal num2str(specimen_num) '_' data_info.ID];

    rm_cmd = ['rm -f ' FigureFileName '.ps'];
    disp(rm_cmd)
    system(rm_cmd)

    disp(['Processing specimen ', num2str(specimen_num)])

    %% plot GIF cleaned map unfiltered
    f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
    t = tiledlayout(length(SN_th_range)/2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
    title(t, 'GIF unfiltered');
    for SN_th = SN_th_range
        nexttile(t);
        data_obj = data_handle_corrected(data_info, data, data_obj.ROI);
        data_obj.generateCleanedDataSamplesGIF(SN_th);
        z = data_obj.read_map(); % Map
        plotElectrodePosition(z,data_obj.ROI,['a= ' num2str(SN_th)],ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    end

    print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
    close(f); % Close the figure to save memory

    %% plot GIF cleaned map filtered
    f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
    t = tiledlayout(length(SN_th_range)/2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
    title(t, 'GIF filtered');
    for SN_th = SN_th_range
        nexttile(t);
        data_obj = data_handle_corrected(data_info, data, data_obj.ROI);
        data_obj.generateCleanedDataSamplesGIF(SN_th);
        z = data_obj.filter_map(data_obj.read_map()); % Map
        plotElectrodePosition(z,data_obj.ROI,['a= ' num2str(SN_th)],ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    end

    print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
    close(f); % Close the figure to save memory


    %% plot GIF JK cleaned map unfiltered
    f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
    t = tiledlayout(length(SN_th_range)/2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
    title(t, 'GIF JK unfiltered');
    for SN_th = SN_th_range
        nexttile(t);
        data_obj = data_handle_corrected(data_info, data, data_obj.ROI);
        data_obj.generateCleanedDataSamplesGIF_JK(SN_th);
        z = data_obj.read_map(); % Map
        plotElectrodePosition(z,data_obj.ROI,['a= ' num2str(SN_th)],ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    end

    print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
    close(f); % Close the figure to save memory

    %% plot GIF cleaned map filtered
    f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
    t = tiledlayout(length(SN_th_range)/2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
    title(t, 'GIF JK filtered');
    for SN_th = SN_th_range
        nexttile(t);
        data_obj = data_handle_corrected(data_info, data, data_obj.ROI);
        data_obj.generateCleanedDataSamplesGIF_JK(SN_th);
        z = data_obj.filter_map(data_obj.read_map()); % Map
        plotElectrodePosition(z,data_obj.ROI,['a= ' num2str(SN_th)],ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    end

    print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
    close(f); % Close the figure to save memory



    disp(' processed and saved to PostScript file.');
    disp(['Output file: ', FigureFileName, '.ps']);

    disp('Convert to pdf')
    cmd = ['ps2pdf ' FigureFileName '.ps ' FigureFileName, '.pdf'];
    disp(cmd)
    system(cmd)

end

function plotElectrodePosition(z,ROI,title_txt,ElectrodePosition,SizeSquare)
    
    %% plot map
    plot_map(z,ROI, 0, 1);
    
    if isempty(ElectrodePosition)
        title(title_txt)
        return
    end

    %% plot electrode position
    size_insert_pix_half= floor(SizeSquare/2);
    
    hold on
    plot(ElectrodePosition(1),ElectrodePosition(2),'x','Color','white')

%     hold on
%     rectangle('Position', [ElectrodePosition(1)-size_insert_pix_half ElectrodePosition(2)-size_insert_pix_half SizeSquare SizeSquare ],'EdgeColor','white')

    %% calculate mean preffered orientation at electrode position
    [X,Y] = meshgrid(1:size(ROI,2),1:size(ROI,1));
    ROI_Electrode = (X >= (round(ElectrodePosition(1))-size_insert_pix_half)).*(X <= (round(ElectrodePosition(1))+size_insert_pix_half)).* ...
                    (Y >= (round(ElectrodePosition(2))-size_insert_pix_half)).*(Y <= (round(ElectrodePosition(2))+size_insert_pix_half));
    
    contour(ROI_Electrode,[1 1],'white','linewidth',1)
    
    z_electrode = z(ROI_Electrode==1);
    meanOrientation = mod(angle(mean(z_electrode,'all'))/2/pi*180,180);%mod(mean(angle(z),'all')/2/pi*180,180);
    title([title_txt ' <ori.> electr.: ' num2str(round(meanOrientation,1)) '°'])

end