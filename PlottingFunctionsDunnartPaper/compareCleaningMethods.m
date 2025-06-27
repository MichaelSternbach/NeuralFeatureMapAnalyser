function compareCleaningMethods(animal,specimen_Num,AnimalDataFolder,FigureFolder,ElectrodePosition)
    %compareCleaningMethods('dunnart',7,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/',[92 110])
    if nargin < 5
        ElectrodePosition = [];
    end
    SizeSquare = 5;
    Bootstrapsamples = 100;
    SN_TH = 4;

    if length(specimen_Num)>1
        %% set Figure file name
        FigureFileName = [FigureFolder 'Cleaningcomparison_' animal];
        rm_cmd = ['rm -f ' FigureFileName '.ps'];
        disp(rm_cmd)
        system(rm_cmd)
    end
    
    for specimen_num = specimen_Num
        %% Load data
        disp('Load data...')
        [data_info, ~, data_obj, data, ~] = getAnimalData(animal, specimen_num, AnimalDataFolder);
        disp('Data loaded successfully.')
        disp('-----------------------')
        
    
        if length(specimen_Num)<2
            %% set Figure file name
            FigureFileName = [FigureFolder 'Cleaningcomparison_' animal num2str(specimen_num) '_' data_info.ID];
            rm_cmd = ['rm -f ' FigureFileName '.ps'];
            disp(rm_cmd)
            system(rm_cmd)
        end
    
        %% set Bootstrapsamples
        data_obj.prepare_samples_array(Bootstrapsamples)
    
        % SNR = calcSNR_OPM_Data(data_obj,false);
    
        % SNR_filtered = calcSNR_OPM_Data(data_obj,true);
    
        %% calc plot limits
        [x_range,y_range] = getRangeXY_ROI(data_obj.ROI);
    
        disp(['Processing specimen ', num2str(specimen_num)])
    
        %% Prepare a new figure for the current specimen
        f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(4,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
    
        title(t,[num2str(specimen_num) '. ' data_info.ID ] )
    
        %% plot not cleaned map unfiltered
        nexttile(t);
        z = data_obj.read_map(); % Map
        
        SNR = calcSNR_OPM_Data(data_obj,false);
        title_txt = ['unfiltered SNR: ' num2str(SNR)];
        plotElectrodePosition(z,data_obj.ROI,title_txt,ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    
        %% plot not cleaned map filtered
        nexttile(t);
        z = data_obj.filter_map(z); % Map
        SNR_filtered = calcSNR_OPM_Data(data_obj,true);
        title_txt = ['filtered SNR: ' num2str(SNR_filtered)];
        plotElectrodePosition(z,data_obj.ROI,title_txt,ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    
        %% plot LSM cleaned map unfiltered
        nexttile(t);
        data_obj.apply_LSM(true);
        z = data_obj.read_map(); % Map
        SNR = calcSNR_OPM_Data(data_obj,false);
        title_txt = ['LSM unfiltered SNR: ' num2str(SNR)];
        plotElectrodePosition(z,data_obj.ROI,title_txt,ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    
        %% plot LSM cleaned map filtered
        nexttile(t);
        z = data_obj.filter_map(z); % Map
        SNR_filtered = calcSNR_OPM_Data(data_obj,true);
        title_txt = ['LSM filtered SNR: ' num2str(SNR_filtered)];
        plotElectrodePosition(z,data_obj.ROI,title_txt,ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    
        %% deactivate LSM
        data_obj.apply_LSM(false);
    
        %% plot GIF cleaned map unfiltered
        nexttile(t);
        data_obj.activateGIF(true,SN_TH)
        z = data_obj.read_map(); % Map
        SNR = calcSNR_OPM_Data(data_obj,false);
        title_txt = ['GIF unfiltered SNR: ' num2str(SNR)];
        plotElectrodePosition(z,data_obj.ROI,title_txt,ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    
        %% plot GIF cleaned map filtered
        nexttile(t);
        z = data_obj.filter_map(z); % Map
        SNR_filtered = calcSNR_OPM_Data(data_obj,true);
        title_txt = ['GIF filtered SNR: ' num2str(SNR_filtered)];
        plotElectrodePosition(z,data_obj.ROI,title_txt,ElectrodePosition,SizeSquare);
        xlim(x_range)
        ylim(y_range)
    
        % %% reset data_obj
        % data_obj = data_handle_corrected(data_info, data, data_obj.ROI);
    
        % %% plot GIF_JK cleaned map unfiltered
        % nexttile(t);
        % data_obj = data_handle_corrected(data_info, data, data_obj.ROI);
        % data_obj.generateCleanedDataSamplesGIF_JK();
        % data_obj.prepare_samples_array(Bootstrapsamples)
        % z = data_obj.read_map(); % Map 
        % SNR = calcSNR_OPM_Data(data_obj,false);
        % title_txt = ['GIF JK unfiltered SNR: ' num2str(SNR)];
        % plotElectrodePosition(z,data_obj.ROI,title_txt,ElectrodePosition,SizeSquare);
        % xlim(x_range)
        % ylim(y_range)
    
        % %% plot GIF_JK cleaned map filtered
        % nexttile(t);
        % z = data_obj.filter_map(z); % Map
        % SNR_filtered = calcSNR_OPM_Data(data_obj,true);
        % title_txt = ['GIF JK filtered SNR: ' num2str(SNR_filtered)];
        % plotElectrodePosition(z,data_obj.ROI,title_txt,ElectrodePosition,SizeSquare);
        % xlim(x_range)
        % ylim(y_range)
    
    
        %% Save the current specimen's plots to a PostScript file
        print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
        close(f); % Close the figure to save memory
    end

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