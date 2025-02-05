function addNewDataSet(mount_point_data,data_info,function_read_data,TestRun)
    if nargin < 4
        TestRun = false;
    end

    if nargin < 2 
        data_info = struct();
    end

    %% data info struct or path to .mat file
    if ~isstruct(data_info)
        if ischar(data_info)
            load(data_info,'data_info');
        else
            error('data_info must be a struct or a path to a .mat file containing a struct')
        end
    end

    %% get data info
    getDataInfo(data_info,mount_point_data);

    %% read data
    if nargin < 3
        if isfile([mount_point_data data_info.preprocessed_data_file])
            load([mount_point_data data_info.preprocessed_data 'data.mat'],'data');
        else
            load(data_info.exp_data_file,'data');

            %% check whether format data is correct and if confirmed overrite data path
            disp(['size data: ' num2str(size(data))])
            disp('data should have the following format: (x-dim y-dim stim block)')
            disp(['stim_order: ' num2str(data_info.stim_order)])
            if isfield(data_info,'field_size_pix')
                disp(['(x-dim y-dim)): ' num2str(data_info.field_size_pix)])
            end
            disp(['exp_data_file: ' data_info.exp_data_file])
            disp(['previous preprocessed_data_file: ' data_info.preprocessed_data_file])
            correctDataformat = input('Is the data format correct and should preprocessed_data_file be set to exp_data_file? (y/n): ','s');
            if correctDataformat
                data_info.preprocessed_data_file = data_info.exp_data_file;
            end
        end
    else
        data = function_read_data(data_info,mount_point_data);

        %% check whether format data is correct and if confirmed save to data path
        disp(['size data: ' num2str(size(data))])
        disp('data should have the following format: (x-dim y-dim stim block)')
        disp(['stim_order: ' num2str(data_info.stim_order)])
        if isfield(data_info,'field_size_pix')
            disp(['(x-dim y-dim)): ' num2str(data_info.field_size_pix)])
        end
        disp(['preprocessed_data_file: ' data_info.preprocessed_data_file])
        if isfile([mount_point_data data_info.preprocessed_data])
            disp('previous data file exists')
        end
        if ~TestRun
            correctDataformat = input('Is the data format correct and should the data be saved to preprocessed_data_file? (y/n): ','s');
            if correctDataformat
                save([mount_point_data data_info.preprocessed_data 'data.mat'],'data','-v7.3');
            end
        end
    end
 

    %% add extra data to data_info
    if ~isfield(data_info,'field_size_pix')
        data_info.field_size_pix = size(data,1:2);
    end

    %% update data_info_file_list (csv)
    updateDataInfoList(data_info,mount_point_data);


    %% set filter settings
    [data_info,power_profile] = setFilterSettingsPowerSpectrum(data_info,data);

    %% make data object
    data_obj = data_handle_corrected(data_info,data);

    %% set ROI
    data_obj = setROI(data_obj);
    ROI = data_obj.ROI;

    if ~TestRun
        %% append ROI to data_file
        save([mount_point_data data_info.preprocessed_data 'data.mat'],'ROI','-append');
    end

    %% fine tune lowpass filter cutoff using the pinwheel numbers
    [data_info,PwPlateus] = fineTuneLowpassFilter(data_info,data_obj,power_profile);

    if ~TestRun

        %% save filter settings
        save([mount_point_data data_info.preprocessed_data 'FilterSettings.mat'],'PwPlateus','power_profile','data_obj')

        %% save figure filter settings
        f = plotPwPlateusAndPowerProfile(PwPlateus,power_profile,partialROIs,data_obj,data_info.settings.lowpass_mm,data_info.settings.highpass_mm,...
                                        [mount_point_data data_info.preprocessed_data 'FilterSettings']);
        close(f)

        %% save data_info to data file
        save([mount_point_data data_info.preprocessed_data 'data.mat'],'data_info','-append');

        %% save data_info as JSON
        writeStructToJSON(data_info,[mount_point_data data_info.preprocessed_data 'data_info.json']);
    end
end

function data_obj = setROI(data_obj)
    %% future features
    % plot CI
    % plot area that encircles areas with CI smaller x
    %apply data cleaning?
    %% use internal function
    data_obj.make_ROI()
end


function [data_info,PwPlateus] = fineTuneLowpassFilter(data_info,data_obj,power_profile)
    %% get OPM
    z = data_obj.filter_map(data_obj.read_map());
    ROI = data_obj.ROI;

    %% initialize figure and parameters
    f = figure;
    partialROIs = [];
    correct_PwPlateu = false;
    lowpass_cutoffs = 0.1:0.1:2;

    disp('For large maps it is recommended to only use one or multiple  smaller ROIs for the pinwheel plateu calculation.')
    while ~correct_PwPlateu

        %% if partialROIs are already set ask if they should be used
        if ~isempty(partialROIs)
            close(f)
            f = plot_partialROIs(z,ROI,partialROIs);
            redoPartialROI = getBooleanInput('Do you want to redo the partial ROIs?');
        else
            redoPartialROI = true;
        end

        if redoPartialROI
            %% plot map and ask for partial ROI
            close(f)
            [partialROIs,f] = getPartialROIs(z,ROI);
        end

        %% check lowpass_cutoffs
        disp('The lowpass_cutoffs to calculate the pinwheel plateus are:') 
        disp([num2str(min(lowpass_cutoffs)) ':'  max(lowpass_cutoffs) 'the mean step size is:' num2str(mean(diff(lowpass_cutoffs)))])
        disp(['the current lowpass_cutoff is: ' num2str(data_info.settings.lowpass_mm)])
        correct_lowpass_cutoffs = getBooleanInput('Do you want to use these lowpass_cutoffs?');
        if ~correct_lowpass_cutoffs
            lowpass_cutoffs = input('Input different lowpass_cutoffs: ');
        end

        %% calculate pinwheel plateus
        PwPlateus = calculatePinwheelPlateus(data_obj,partialROIs);
        close(f)
        f = plotPwPlateus(PwPlateus,data_info.settings.lowpass_mm);
        correct_PwPlateu = getBooleanInput('Are the pinwheel plateus usable or do you want to redo them ?');
    end

    %% fine tune lowpass filter cutoff
    disp('Fine tune lowpass filter cutoff using the pinwheel numbers')
    disp('Lowpass filter cutoff (mm) should be inside the pinwheel plateus')
    correctFilter = false;
    while ~correctFilter
        lowpass_mm = getNumericInput('lowpass_mm: ');
        close(f)
        f = plotPwPlateusAndPowerProfile(PwPlateus,power_profile,partialROIs,data_obj,lowpass_mm,data_info.settings.highpass_mm);
        correctFilter = getBooleanInput('Is the lowpass cutoff correct?');
    end
    close(f)

    data_info.settings.lowpass_mm = lowpass_mm;
    
end

function f = plotPwPlateusAndPowerProfile(PwPlateus,power_profile,partialROIs,data_obj,lowpass_mm,highpass_mm,SavFig)
    if nargin < 7
        SavFig = '';
    end


    f = figure;
    t = tiledlayout(3,1);

    %% plot power profile
    nexttile(t)
    plot(power_profile.scale_mm,power_profile.values);
    hold on
    plot([lowpass_mm lowpass_mm],[0 max(power_profile.values)],'r','DisplayName','lowpass cutoff')
    hold on
    plot([highpass_mm highpass_mm],[0 max(power_profile.values)],'r','DisplayName','highpass cutoff')
    xlabel('Scale in mm')
    ylabel('Power')
    legend()
    xlim([0 max(profile_range_mm)])
    title('Power profile')

    %% plot pinwheel plateus
    color_list = getColorList(length(PwPlateus));
    nexttile(t)
    for ii = 1:length(PwPlateus)
        plot(PwPlateus{ii}.global_plateau.lowpass_vs_density(:,1),PwPlateus{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)],'Color',color_list{ii})
        hold on
    end
    plot([lowpass_mm lowpass_mm],[0 max(PwPlateus.values)],'r','DisplayName','current lowpass cutoff')

    %% plot ROIs
    nexttile(t)
    z = data_obj.filter_map(data_obj.read_map());
    [x_range,y_range] = getRangeXY_ROI(ROI);
    ROI = data_obj.ROI;
    plot_map(z,ROI,0,1)
    for ii = 1:size(partialROIs,3)
        contour(partialROIs(:,:,ii),[1 1],'LineWidth',2,'Color',color_list{ii});
    end
    xlim(x_range)
    ylim(y_range)
    title('ROIs')

    if ~isempty(SavFig)
        print(f, '-dpsc', '-fillpage', [SaveFig, '.ps']);
        close(f);
    end
end

function color_list = getColorList(n)
    color_list = cell(n,1);
    for ii = 1:n
        color_list{ii} = [rand rand rand];
        for jj = 1:ii-1
            while sum(abs(color_list{ii}-color_list{jj})) < 0.5
                color_list{ii} = [rand rand rand];
            end
        end
    end
end 

function f = plotPwPlateus(PwPlateus,lowpass_cutoff)
    f = figure;
    for ii = 1:length(PwPlateus)
        plot(PwPlateus{ii}.global_plateau.lowpass_vs_density(:,1),PwPlateus{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)])
        hold on
    end
    plot([lowpass_cutoff lowpass_cutoff],[0 max(PwPlateus.values)],'r','DisplayName','current lowpass cutoff')
    xlabel('lowpass cutoff in mm')
    ylabel('density')
    legend()
end

function PwPlateus = calculatePinwheelPlateus(data_obj,partialROIs,lowpass_cutoffs)

    PwPlateus = cell(length(splitROI),1);
    for ii = 1:length(splitROI)
        data_obj.set_ROI(partialROIs(:,:,ii));
        PwPlateus{ii} = find_lowpassPwNumber(data_obj,data_info,highpass_mm,lowpass_mm,lowpass_cutoffs);
        % plot(PwPlateus{ii}.global_plateau.lowpass_vs_density(:,1),PwPlateus{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)],'Color',color_list{ii})
        % hold on
    end

end


function [partialROIs,f] = getPartialROIs(z,ROI)
    f = figure;
    plot_map(z,ROI,0,1)
    doPartialROI = getBooleanInput('Do you want to use a partial ROI? ');
    if doPartialROI
        disp('Select partial ROIs by hand using the format [[x1 y1 x2 y2];[x1 y1 x2 y2]] or input a number of ROIs to select automatically!')

        CorrectPartialROIs = false;
        while ~CorrectPartialROIs
            %% convert partialROIs to uniform format
            SplitROI = getNumericInput('Input partialROIs: ');
            partialROIs = calculatePartialROIs(SplitROI);

            close(f)
            f = plot_partialROIs(z,ROI,partialROIs);
            CorrectPartialROIs = getBooleanInput('Do you want to use these partial ROIs? ');
        end
    else
        partialROIs = reshape(ROI,[size(ROI,1),size(ROI,2),1]);
    end
end

function f = plot_partialROIs(z,ROI,partialROIs)
    [x_range,y_range] = getRangeXY_ROI(ROI);
    f = figure;
    plot_map(z,ROI,0,1)
    hold on
    for ii = 1:size(partialROIs,3)
        contour(partialROIs(:,:,ii),[1 1],'white','LineWidth',2);
    end
    xlim(x_range)
    ylim(y_range)
end

function partialROIs = calculatePartialROIs(SplitROI)
    if length(splitROI)>1
        partialROIs = calcPartialROIs(data_obj.ROI,splitROI);
    else
        partialROIs = calcSplitROI(data_obj.ROI,SplitROI);
    end

end

function partialROIs = calcSplitROI(ROI,NumSplitROI)
    %% get ranges ROI
    [x_range,y_range] = getRangeXY_ROI(ROI);
    %% get split ROI
    x_split = round(linspace(x_range(1),x_range(2),NumSplitROI+1));
    y_split = round(linspace(y_range(1),y_range(2),NumSplitROI+1));
    size_splitROI = (x_split(2)-x_split(1))*(y_split(2)-y_split(1));
    n=1;
    for ii = 1:NumSplitROI
        for jj = 1:NumSplitROI
            %% get splitROI
            splitROI_x = [x_split(ii) x_split(ii+1)];
            splitROI_y = [y_split(jj) y_split(jj+1)];

            %% calc ROI
            partialROI = calcPartialROI(splitROI_x,splitROI_y,ROI);

            %% check size of splitROI
            if (sum(partialROI,'all')/size_splitROI) > 0.8            
                % splitROI{n}.x = splitROI_x;
                % splitROI{n}.y = splitROI_y;
                partialROIs(:,:,n) = partialROI;
                n = n+1;
            end
        end
    end

end

function partROI = calcPartialROI(x,y,ROI)
    %% get Meshrid size ROI
    [X,Y] = meshgrid(1:size(ROI,2),1:size(ROI,1));
    partROI_bool = ((X >= x(1) & X <= x(2)) & (Y >= y(1) & Y <= y(2)));
    partROI = zeros(size(ROI));
    partROI(partROI_bool)=1;
    partROI = (partROI.*ROI);
end                                                                                                                                                                                                                

function partialROIs = calcPartialROIs(ROI,splitROI)
    partialROIs = zeros(size(ROI,1),size(ROI,2),length(splitROI));
    for ii = 1:size(splitROI,1)
        partialROIs(:,:,ii) = calcPartialROI(splitROI(ii,[1 3]),splitROI(ii,[2 4]),ROI);
    end
end

%contour(A, [0.5 0.5], 'r', 'LineWidth', 2);

function [data_info,power_profile] = setFilterSettingsPowerSpectrum(data_info,data)

    %% parameters for profile
    max_dim_mm = max(data_info.field_size_pix)/data_info.pix_per_mm;
    profile_range_mm = 0.1:0.01:0.8*max_dim_mm;
    
    %% calculate power profile

    GoodPowerprofile = false;
    f= figure();
    while ~GoodPowerprofile
        close(f)
        power_profile = define_filter_settings(data_info,ones(size(data,1:2)),data_obj.data,profile_range_mm);
        f = plotPowerProfile(power_profile);
        GoodPowerprofile = getBooleanInput('Do you want to use this Powerprofile?');
        if ~GoodPowerprofile
            profile_range_mm = input('Input different profile_range_mm (0.1:0.01:5): ');
        end
    end

    disp('Set filter settings based on power profile!')
    disp('Lowpass filter cutoff (mm) should be smaller than the peak of the main signal')
    disp('Highpass filter cutoff (mm) should be larger than the peak of the main signal')
    lowpass_mm = getNumericInput('lowpass_mm: ');
    highpass_mm = getNumericInput('highpass_mm: ');
    hold on
    plot([lowpass_mm lowpass_mm],[0 max(power_profile.values)],'r','DisplayName','lowpass cutoff')
    hold on
    plot([highpass_mm highpass_mm],[0 max(power_profile.values)],'r','DisplayName','highpass cutoff')
    legend()

    saveFilter = false;
    while ~saveFilter
        saveFilter = input('Save filter settings? (y/n): ','s');
        if saveFilter
            data_info.settings.lowpass_mm = lowpass_mm;
            data_info.settings.highpass_mm = highpass_mm;
            close(f)
        else
            lowpass_mm = getNumericInput('lowpass_mm: ');
            highpass_mm = getNumericInput('highpass_mm: ');
            close(f)
            f = plotPowerProfile(power_profile);
            hold on
            plot([lowpass_mm lowpass_mm],[0 max(power_profile.values)],'r','DisplayName','lowpass cutoff')
            hold on
            plot([highpass_mm highpass_mm],[0 max(power_profile.values)],'r','DisplayName','highpass cutoff')
            legend()
        end
    end
end

function f = plotPowerProfile(power_profile)
    f = figure;
    plot(power_profile.scale_mm,power_profile.values);
    xlabel('Scale in mm')
    ylabel('Power')
    xlim([0 max(profile_range_mm)])
    set(gca,'fontsize',15)
end



function updateDataInfoList(data_info)
    %% get dir of this script
    script_dir = fileparts(mfilename('fullpath'));

    %% pull latest version of script_dir from git
    system(['git pull ' script_dir]);

    %% get data_info_file_list path
    data_info_file_list_file_path = [script_dir data_info.data_info_file_list];

    %% read in data_info_file_list
    data_info_file_list = readtable(data_info_file_list_file_path);

    %% get values Dataset,MakeInfo,DestinationFolder,ExperimentID,ExperimentName
    Dataset = data_info.animal;
    MakeInfo = data_info.data_info_file_path;
    
    if isfield(data_info,'DestinationFolder')
        DestinationFolder = data_info.DestinationFolder;
    else
        DestinationFolder = '';
    end

    ExperimentID = data_info.exp_num;
    ExperimentName = data_info.ID

    %% check whether dataset already exists
    if any(strcmp(data_info_file_list.Dataset,Dataset) & strcmp(data_info_file_list.ExperimentID,ExperimentID))
        disp(['Dataset ' Dataset ' with ExperimentID ' num2str(ExperimentID) ' already exists in data_info_file_list'])
        %% disp full MakeInfo and ExperimentName in table
        disp(data_info_file_list(strcmp(data_info_file_list.Dataset,Dataset) & strcmp(data_info_file_list.ExperimentID,ExperimentID),:))

    else
        disp(['Adding new dataset ' Dataset ' with ExperimentID ' num2str(ExperimentID) ' to data_info_file_list'])
        %% add new row to table. columns: Dataset,MakeInfo,DestinationFolder,ExperimentID,ExperimentName
        data_info_file_list = [data_info_file_list;{Dataset,MakeInfo,DestinationFolder,ExperimentID,ExperimentName}];

        %% write table to data_info_file_list
        writetable(data_info_file_list,data_info_file_list_file_path);

        %% commit changes to git
        system(['git add ' data_info_file_list_file_path ; 'git commit -m "added new dataset enty ' Dataset num2str(ExperimentID) '"']);
        disp('Changes saved and commited to git')
    end  

end