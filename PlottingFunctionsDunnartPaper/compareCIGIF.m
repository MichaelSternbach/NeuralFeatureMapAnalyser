function compareCIGIF(animal,specimen_num_list,AnimalDataFolder,ResultDataFolder,FigureFolder,DoFilter,Bootstrapsamples)
    %close all; compareCIGIF('dunnart',6,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/DataHPC/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/')
    if nargin <6
        DoFilter = false;
    end
    if nargin <7
        Bootstrapsamples = 100;
    end

    alpha = 0.05;


    if length(specimen_num_list) >1
        %% set Figure file name
        FigureFileName = [FigureFolder 'ComparisonCIGIF_' animal '_all'];
        rm_cmd = ['rm -f ' FigureFileName '.ps'];
        disp(rm_cmd)
        system(rm_cmd)
    end
    
    for specimen_num = specimen_num_list

        disp(['Processing specimen ', num2str(specimen_num)])


        %% Load data
        disp('Load data...')
        [data_info, data_path, data_obj, data, BloodVesselImg] = getAnimalData(animal, specimen_num, AnimalDataFolder);
        disp('Data loaded successfully.')
        disp('-----------------------')

        data_obj.prepare_samples_array(Bootstrapsamples)

        %% calc plot limits
        [x_range,y_range] = getRangeXY_ROI(data_obj.ROI);

        if length(specimen_num_list) == 1
            %% set Figure file name
            FigureFileName = [FigureFolder 'ComparisonCIGIF_' animal num2str(specimen_num) '_' data_info.ID];

            rm_cmd = ['rm -f ' FigureFileName '.ps'];
            disp(rm_cmd)
            system(rm_cmd)
        end

        %% plot GIF cleaned map unfiltered
        f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
        title(t, ['CI Comparison ' data_info.ID])

        ax = nexttile(t);
        DataCleaning = 'none';
        data_obj.activateGIF(false)
        data_obj.apply_LSM(false);
        %data_obj.generateCleanedDataSamplesGIF();
        DataFolder = [ResultDataFolder animal '/' animal num2str(specimen_num) '/' DataCleaning '/' num2str(Bootstrapsamples) '/'];
        mkdir(DataFolder)
        [CI_angle,CI_Abs,ROI] = getCI(data_obj,alpha,'bca',DoFilter);
        plot_mapAbs(CI_angle,'',180,0,data_obj.ROI,ax)
        title('None')
        xlim(x_range)
        ylim(y_range)


        ax = nexttile(t);
        DataCleaning = 'LSM';
        data_obj.activateGIF(false)
        data_obj.apply_LSM(true);
        DataFolder = [ResultDataFolder animal '/' animal num2str(specimen_num) '/' DataCleaning '/' num2str(Bootstrapsamples) '/'];
        mkdir(DataFolder)
        CI = calcCIs(data_obj,alpha,DoFilter,DataFolder);
        CI_angle = CI.BCA.CI_angle;
        plot_mapAbs(CI_angle,'',180,0,data_obj.ROI,ax)
        title('LSM')
        xlim(x_range)
        ylim(y_range)

        ax = nexttile(t);
        DataCleaning = 'GIF';
        data_obj.activateGIF(true,3)
        data_obj.apply_LSM(false);
        DataFolder = [ResultDataFolder animal '/' animal num2str(specimen_num) '/' DataCleaning '/' num2str(Bootstrapsamples) '/'];
        mkdir(DataFolder)
        [CI.BCA.CI_angle,CI.BCA.CI_Abs,CI.BCA.ROI] = getCI(data_obj,alpha,'bca',DoFilter);
        CI_angle = CI.BCA.CI_angle;
        plot_mapAbs(CI_angle,'',180,0,data_obj.ROI,ax)
        title('GIF')
        xlim(x_range)
        ylim(y_range)


%         ax = nexttile(t);
%         DataCleaning = 'LSM+GIF';
%         data_obj.activateGIF(true,3)
%         data_obj.apply_LSM(true);
%         DataFolder = [ResultDataFolder animal '/' animal num2str(specimen_num) '/' DataCleaning '/' num2str(Bootstrapsamples) '/'];
%         mkdir(DataFolder)
%         [CI.BCA.CI_angle,CI.BCA.CI_Abs,CI.BCA.ROI] = getCI(data_obj,alpha,'bca',DoFilter);
%         CI_angle = CI.BCA.CI_angle;
%         plot_mapAbs(CI_angle,'',180,0,data_obj.ROI,ax)
%         title('LSM+GIF')
%         xlim(x_range)
%         ylim(y_range)

        print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
        close(f); % Close the figure to save memory

    end
%     ax = nexttile(t);
%     DataCleaning = 'GIF_JK';
%     data_obj = data_handle_corrected(data_info, data, data_obj.ROI);
%     data_obj.generateCleanedDataSamplesGIF_JK();
%     data_obj.prepare_samples_array(Bootstrapsamples)
%     DataFolder = [ResultDataFolder animal '/' animal num2str(specimen_num) '/' DataCleaning '/' num2str(Bootstrapsamples) '/'];
%     mkdir(DataFolder)
%     CI = calcCIs(data_obj,alpha,DoFilter,DataFolder);
%     CI_angle = CI.BCA.CI_angle.*sqrt(size(data_obj.data,4));
%     CI_angle(find(CI_angle>180))=180;
%     plot_mapAbs(CI_angle,'',180,0,data_obj.ROI,ax)
%     title('GIF JK*sqrt(#samples)')
%     xlim(x_range)
%     ylim(y_range)



    



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