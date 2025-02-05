function testCleaningSteps(animal,specimen_Num,AnimalDataFolder,FigureFolder,ElectrodePosition)
    %testCleaningSteps('dunnart',7,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/',[92 110])
    if nargin < 5
        ElectrodePosition = [];
    end
    SizeSquare = 5;
    Bootstrapsamples = 100;
    profile_range_mm = 0.2:0.01:2.5;
    SN_TH=4;
    alpha = 0.05;

    if length(specimen_Num)>1
        %% set Figure file name
        FigureFileName = [FigureFolder 'CleaningSteps_' animal];
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
        

        %% set filter settings
        lowpass_mm = 0.4;
        highpass_mm = 2;

        data_info.settings.lowpass_mm = lowpass_mm;
        data_info.settings.highpass_mm = highpass_mm;

        data_obj.set_filter_parameters("lowpass",lowpass_mm);
        data_obj.set_filter_parameters("highpass",highpass_mm);
    
        if length(specimen_Num)<2
            %% set Figure file name
            FigureFileName = [FigureFolder 'CleaningSteps_' animal num2str(specimen_num) '_' data_info.ID];
            rm_cmd = ['rm -f ' FigureFileName '.ps'];
            disp(rm_cmd)
            system(rm_cmd)
        end
    
        %% calculate SNR
        data_obj.prepare_samples_array(Bootstrapsamples)
    
        % SNR = calcSNR_OPM_Data(data_obj,false);
    
        % SNR_filtered = calcSNR_OPM_Data(data_obj,true);
    
        %% calc plot limits    
        disp(['Processing specimen ', num2str(specimen_num)])
    
        %% Prepare a new figure for the current specimen
        f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(4,3, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
    
        title(t,[num2str(specimen_num) '. ' data_info.ID ] )
    
        %% plot not cleaned map unfiltered
        z = data_obj.read_map(); % Map
        %NoiseMeasure.SNR = calcSNR_OPM_Data(data_obj,false);
        [NoiseMeasure.CI_angle,~,~] = getCI(data_obj,alpha,'bca',false);
        plotMapSignal(t,z,data_obj.ROI,ElectrodePosition,SizeSquare,data_info,NoiseMeasure,profile_range_mm,"Not cleaned unfiltered")
        

        %% plot not cleaned map filtered
        z = data_obj.filter_map(z); % Map
        %NoiseMeasure.SNR = calcSNR_OPM_Data(data_obj,true);
        [NoiseMeasure.CI_angle,~,~] = getCI(data_obj,alpha,'bca',true);
        plotMapSignal(t,z,data_obj.ROI,ElectrodePosition,SizeSquare,data_info,NoiseMeasure,profile_range_mm,"Not cleaned filtered")


        % %% plot LSM cleaned map unfiltered
        % data_obj.apply_LSM(true);
        % z = data_obj.read_map(); % Map
        % SNR = calcSNR_OPM_Data(data_obj,false);
        % plotMapSignal(t,z,data_obj.ROI,ElectrodePosition,SizeSquare,data_info,SNR,profile_range_mm,"LSM cleaned unfiltered")
        % data_obj.apply_LSM(false);

        %% plot GIF cleaned map unfiltered
        data_obj.activateGIF(true,SN_TH)
        z = data_obj.read_map(); % Map
        %NoiseMeasure.SNR = calcSNR_OPM_Data(data_obj,false);
        [NoiseMeasure.CI_angle,~,~] = getCI(data_obj,alpha,'bca',false);
        plotMapSignal(t,z,data_obj.ROI,ElectrodePosition,SizeSquare,data_info,NoiseMeasure,profile_range_mm,"GIF cleaned unfiltered")

        %% plot GIF cleaned map filtered
        z = data_obj.filter_map(z); % Map
        %NoiseMeasure.SNR = calcSNR_OPM_Data(data_obj,true);
        [NoiseMeasure.CI_angle,~,~] = getCI(data_obj,alpha,'bca',true);
        plotMapSignal(t,z,data_obj.ROI,ElectrodePosition,SizeSquare,data_info,NoiseMeasure,profile_range_mm,"GIF cleaned filtered")
    
    
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

function plotMapSignal(t,z,ROI,ElectrodePosition,SizeSquare,data_info,NoiseMeasure,profile_range_mm,title_txt)

    %% plot map
    nexttile(t);
    [x_range,y_range] = getRangeXY_ROI(ROI);
    plotElectrodePosition(z,ROI,'Signal',ElectrodePosition,SizeSquare)
    xlim(x_range)
    ylim(y_range)
    title(title_txt)

    %% plot powerspetrum
    nexttile(t);
    power_profile = define_filter_settings(data_info,ROI,z,profile_range_mm);
    plot(power_profile.k_mm_inv,power_profile.values_kspace)
    hold on
    % plot filter parameters
    plot([1/data_info.settings.lowpass_mm 1/data_info.settings.lowpass_mm],[0 max(power_profile.values_kspace)],"DisplayName","Lowpass")
    plot([1/data_info.settings.highpass_mm 1/data_info.settings.highpass_mm],[0 max(power_profile.values_kspace)],"DisplayName","Highpass")
    xlabel('k [1/mm]')
    ylabel('Power')
    title('  Power Spectrum')
    %legend()
    
    %% plot CI
    ax = nexttile(t);
    plot_mapAbs(NoiseMeasure.CI_angle,'CI angle [°]',180,0,ROI,ax)
    xlim(x_range)
    ylim(y_range)

%     %% plot text signal and noise power + SNR
%     nexttile(t);
%     SNR = NoiseMeasure.SNR
%     SignalPower = mean(z(ROI).*conj(z(ROI)),'all');
%     NoisePower = SignalPower/SNR;
%     text(0.1,0.9,['Signal Power: ' num2str(SignalPower)],'Units','normalized')
%     text(0.1,0.8,['Noise Power: ' num2str(NoisePower)],'Units','normalized')
%     text(0.1,0.7,['SNR: ' num2str(SNR)],'Units','normalized')
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
    

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