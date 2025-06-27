function compareParameterGIF_DirectionData(animal,specimen_num_list,AnimalDataFolder,FigureFolder,ElectrodePosition,SN_th_range)
    %compareParameterGIF('dunnart',7,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/',[92 110])
    %close all; compareParameterGIF('dunnart',6,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/',[67 96])
    if nargin < 5
        ElectrodePosition = [];
    end
    if nargin <6
        %SN_th_range = 0.5:0.5:4;
        SN_th_range = [0 0.001 0.01 0.1 3 10 100 1000];
    end

    SN_th_range_long = [ 0 0.0001 0.001 0.01 0.1 3 10 100 1000 5000];%[0 0.000001 0.00001 0.0001 0.001 0.01 0.1 1 3 10 20 50 100 500 1000 2000 5000 ];

    SizeSquare = 5;
    Bootstrapsamples = 100;
    alpha = 0.05;

    if length(specimen_num_list)>1
        %% set Figure file name
        FigureFileName = [FigureFolder 'ComparisonParameterGIF_' animal '_DirectionData'];
    
        %% remove old file
        rm_cmd = ['rm -f ' FigureFileName '.ps'];
        disp(rm_cmd)
        system(rm_cmd)
    end

    %% loop over specimens
    for specimen_num = specimen_num_list
        disp(['Processing specimen ', num2str(specimen_num)])

        %% Load data
        disp('Load data...')
        [data_info, data_path, data_obj, data, BloodVesselImg] = getAnimalData(animal, specimen_num, AnimalDataFolder);

        [~,~,DirectionData] = getDirectionData(data_info,data_path,0);
        DirectionData  = real(DirectionData);

        %% make data obj dir data
        stimDir = [data_info.stim_order(find(~isnan(data_info.stim_order))) data_info.stim_order(find(~isnan(data_info.stim_order)))+180];
        stimDir = [stimDir NaN];
        data_info.stim_order=stimDir;
        data_obj = data_handle_corrected(data_info,DirectionData,data_obj.ROI);
        data_obj.prepare_samples_array(Bootstrapsamples)

        disp('Data loaded and prepared successfully.')
        disp('-----------------------')

        if length(specimen_num_list)==1
            %% set Figure file name
            FigureFileName = [FigureFolder 'ComparisonParameterGIF_' animal num2str(specimen_num) '_' data_info.ID 'DirectionData'];

            %% remove old file
            rm_cmd = ['rm -f ' FigureFileName '.ps'];
            disp(rm_cmd)
            system(rm_cmd)
        end

        %% calc plot limits
        [x_range,y_range] = getRangeXY_ROI(data_obj.ROI);

        %% get SNR and dimensions of 
        %[data_clean,phi,SN,sn] = generateCleanedDataSamplesGIF(obj,SN_th)

        SN_list = zeros(length(SN_th_range_long),1);
        sn_max_list = zeros(length(SN_th_range_long),1);
        sn_size_list = zeros(length(SN_th_range_long),1);

        for ii = 1:length(SN_th_range_long)
            data_obj.activateGIF(true,SN_th_range_long(ii))
            [~, out] = data_obj.read_map();

    %         SNR = calcSNR_OPM_Data(data_obj,false);
    %         SN_list(i) = SNR;
            SN_list(ii) = mean(out.SN.*squeeze(mean(abs(out.phi),[1 2])))/mean(abs(out.phi),'all');
            sn_max_list(ii) = max(real(diag(out.sn)));
            sn_size_list(ii) = sum(diag(out.sn)>0,"all");
        end

        f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
        title(t, ['GIF parameters ' data_info.ID]);

        %% plot SNR
        nexttile(t);
        loglog(SN_th_range_long,SN_list)
        title('SNR')
        xlabel('SN threshold')
        ylabel('weighted Mean SNR GIFs')

        %% plot max sn
        nexttile(t);
        semilogx(SN_th_range_long,sn_max_list)  
        title('max EV')
        xlabel('SN threshold')
        ylabel('max EV')

        %% plot size solution space
        nexttile(t);
        semilogx(SN_th_range_long,sn_size_list)
        title('size solution space')
        xlabel('SN threshold')
        ylabel('dim solution space')

        print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
        close(f); % Close the figure to save memory


        %% plot GIF cleaned map unfiltered
        f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(length(SN_th_range)/2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
        title(t, ['GIF unfiltered ' data_info.ID]);
        for SN_th = SN_th_range
            nexttile(t);
            data_obj.activateGIF(true,SN_th)
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
        title(t, ['GIF filtered ' data_info.ID]);
        for SN_th = SN_th_range
            nexttile(t);
            data_obj.activateGIF(true,SN_th)
            z = data_obj.filter_map(data_obj.read_map()); % Map
            plotElectrodePosition(z,data_obj.ROI,['a= ' num2str(SN_th)],ElectrodePosition,SizeSquare);
            xlim(x_range)
            ylim(y_range)
        end
    
        print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
        close(f); % Close the figure to save memory

        %% plot CI of GIF cleaned map unfiltered
        f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(length(SN_th_range)/2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
        title(t, ['CI GIF unfiltered ' data_info.ID]);
        for SN_th = SN_th_range
            ax=nexttile(t);
            data_obj.activateGIF(true,SN_th)
            DoFilter = false;
            [CI_angle,CI_Abs,ROI] = getCI(data_obj,alpha,'bca',DoFilter);
            plot_mapAbs(CI_angle,['a= ' num2str(SN_th)],180,0,data_obj.ROI,ax)
            xlim(x_range)
            ylim(y_range)
        end

        print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
        close(f); % Close the figure to save memory

        %% plot CI of GIF cleaned map filtered
        f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(length(SN_th_range)/2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
        title(t, ['CI GIF filtered ' data_info.ID]);
        for SN_th = SN_th_range
            ax=nexttile(t);
            data_obj.activateGIF(true,SN_th)
            DoFilter = true;
            [CI_angle,CI_Abs,ROI] = getCI(data_obj,alpha,'bca',DoFilter);
            plot_mapAbs(CI_angle,['a= ' num2str(SN_th)],180,0,data_obj.ROI,ax)
            xlim(x_range)
            ylim(y_range)
        end

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