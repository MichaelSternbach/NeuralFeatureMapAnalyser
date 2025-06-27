function [data_info,data_obj,data] = ProcessDataOPM(data,data_info,ROI,result_dir,removeNanStimSignal,number_bootstrapsamples)
    %% ProcessDataOPM
    % This function processes OPM data and saves the results to a specified location.
    % data: OPM data to be processed (x-dim, y-dim, stim, block) or (x-dim, y-dim, stim, block, time).
    % data_info: struct containing information about the data. 
    % Minimum fields: stim_order, pix_per_mm, animal, ID.
    % if #-dim data is 5, the data_info struct must contain the field sigWin.
    % Optional (but encouraged) fields: weight_in_grams, age_days, gender, date_recording
    % ROI: region of interest to be used for processing. size: (x-dim, y-dim)
    % In case data, data_info and ROI are provided as strings, the function will load them from the specified path.
    % result_dir: directory where the results will be saved.
    % If not provided, the function will create a directory in the current working directory with the name of the animal and ID.
    % removeNanStimSignal: boolean indicating whether to remove stimulus independent signal.
    % number_bootstrapsamples: number of bootstrap samples to be used for the analysis.

    %% add path to OPM_Pocessing functions
    addpath 'OPM_Processing' 
    
    %% set default values
    if nargin < 2 
        data_info = struct();
    end
    if (nargin < 3 && isempty(ROI) && ~isfield(data_info,'ROI'))
        disp('No ROI provided, using all data as ROI')
        ROI = ones(size(data,1),size(data,2));
        % add make ROI option
    end
    if nargin <5
        removeNanStimSignal = false;
    end
    if nargin < 6
        number_bootstrapsamples = 100;
    end

    %% check format data and update data_info
    [data_info,data] = checkDataFormats(data_info,data,ROI);

    %% define and make result_dir if it does not exist
    if nargin < 4
        result_dir = makeResultDir(data_info);
    end
    if ~isfolder(result_dir)
        mkdir(result_dir);
    end


    
    if removeNanStimSignal
        %% compare methods for the removal of stimulus independent signal
        displayHeader('Compare methods for the removal of stimulus independent signal')
        [data_info,data] = compareMethodsRemovalStimulusIndependentSignal(data_info,data,ROI,result_dir,number_bootstrapsamples);
    else
        %% if data is 5D, average over sigWin
        if length(size(data)) == 5
            data = mean(data(:,:,:,:,data_info.sigWin),5);
        end

        %% save data as .mat
        save([result_dir 'ProcessedData.mat'],'data');
    end

    

    %% prepare data object
    displayHeader('Prepare data object')
    data_obj = prepareDataObject(data_info,data,ROI,number_bootstrapsamples);

    %% if no ROI was not provided, set it based on signal
    if any(ROI== ones(size(ROI)))
        displayHeader('Determine ROI from data')
        data_obj = determineROIFromData(data_obj);
    end

    %% compare data cleaning methods
    displayHeader('Compare data cleaning methods')
    [data_info,data_obj] = compareDataCleaningMethods(data_info,data_obj,result_dir);
    

    %% set filter settings based on spacial analysis results
    displayHeader('Set filter settings based on spacial analysis results')
    [data_info,data_obj,spacial_analysis_results] = setFilterSettingsBasedOnSpacialAnalysis(data_info,data_obj,result_dir);

    %% do pinwheel density analysis to refine lowpass filter cutoff
    displayHeader('Refine lowpass filter cutoff based on pinwheel density analysis')
    [data_info,data_obj] = refineLowpassFilterBasedOnPinwheelDensityAnalysis(data_info,data_obj,spacial_analysis_results,result_dir);


    %% save data info as .mat and .json 
    displayHeader('Save data info')
    ROI = data_obj.ROI;
    save([result_dir 'data_info.mat'],'data_info','data_obj','ROI');
    writeStructToJSON(data_info,[result_dir 'data_info.json']);

    disp('Processing finished!')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare function variables
function [data_info,data]= checkDataFormats(data_info,data,ROI)
    %% data info struct or path to .mat file
    if ~isstruct(data_info)
        if ischar(data_info)
            load(data_info,'data_info');
        else
            error('data_info must be a struct or a path to a .mat file containing a struct')
        end
    end

    %% data array or path to .mat file
    if ~isnumeric(data)
        if ischar(data)
            load(data,'data');
        else
            error('data must be a numeric array or a path to a .mat file containing a numeric array')
        end
    end
    %% ROI array or path to .mat file
    if ~isnumeric(ROI) & ~islogical(ROI)
        if ischar(ROI)
            load(ROI,'ROI');
        else
            error('ROI must be a numeric array or a path to a .mat file containing a numeric array')
        end
    end

    %% check format data info 
    if ~isfield(data_info,'stim_order')
        disp('data_info must contain stim_order field!')
        stim_order = getNumericInput('Input stim_order: ');
        data_info.stim_order = stim_order;
    end
    if ~isfield(data_info,'pix_per_mm')
        disp('data_info must contain pix_per_mm field!')
        pix_per_mm = getNumericInput('Input pix_per_mm: ');
        data_info.pix_per_mm = pix_per_mm;
    end 
    if ~isfield(data_info,'animal')
        disp('data_info must contain animal field!')
        animal = getStringInput('Input animal: ');
        data_info.animal = animal;
    end
    if ~isfield(data_info,'ID')
        disp('data_info must contain ID field!')
        ID = getStringInput('Input ID: ');
        data_info.ID = ID;
    end

    %% check format data
    if length(size(data)) < 4 || length(size(data)) > 5
        error('data must be 4D or 5D')
    end
    if length(size(data)) == 5 & ~isfield(data_info,'sigWin')
        error('data_info must contain sigWin field if data is 5D')
    end
    if size(data,3) ~= length(data_info.stim_order)
        error('data must have the same number of stim as data_info.stim_order')
    end

    %% check format ROI and add to data_info
    if ~isfield(data_info,'ROI')
        data_info.ROI = ROI;
    end
    if size(data_info.ROI,1) ~= size(data,1) || size(data_info.ROI,2) ~= size(data,2)
        error('ROI must have the same size as data')
    end
    if size(data_info.ROI,3) ~= 1
        error('ROI must be a 2D matrix')
    end


    %% add field size to data_info
    if ~isfield(data_info,'field_size_pix')
        data_info.field_size_pix = size(data,1:2);
    end
    
end

function result_dir = makeResultDir(data_info)
    executing_dir = fileparts(mfilename('fullpath'));
    main_result_dir = [executing_dir '/' data_info.animal '/'];
    if ~isfolder(main_result_dir)
        mkdir(main_result_dir);
    end
    result_dir = [main_result_dir data_info.ID '/'];
end

%%% main processing steps
function displayHeader(text)

    %% make headertext
    header_txt = ['## ' text ' ##'];
    
    %% make row of # same length as header_txt
    row_txt = repmat('#',1,length(header_txt));

    %% display header
    disp(row_txt)
    disp(header_txt)
    disp(row_txt)
end

function data_obj = determineROIFromData(data_obj)
    %% set ROI based on data
    data_obj.make_ROI();
    % %% extract ROI
    % ROI = data_obj.ROI;
end

function [data_info,data] = compareMethodsRemovalStimulusIndependentSignal(data_info,data,ROI,result_dir,number_bootstrapsamples)
    displayHeader('Compare methods for the removal of stimulus independent signal')
    
    result_file = [result_dir 'nan_stim_removal_methods_data.mat'];

    if exist(result_file,'file')
        %% load data
        load(result_file,'data_removal_methods');
        disp('Data for the removal of stimulus independent signal already exist, loading them from file')
        disp(result_file)
    else
        %% set removal_method_list
        if length(size(data)) == 5
            removal_method_list = {'none','ff','mean','cocktail party'};
        elseif length(size(data)) == 4
            removal_method_list = {'none','cocktail party','mean'};
        else
            error('Data must be 4D or 5D')
        end
        %% calculate results for different methods
        data_removal_methods = calcRemovalMethodsComparisonData(data,ROI,removal_method_list,number_bootstrapsamples);
        
        %% save data
        save(result_file,'data_removal_methods');
    end

    %% plot figures
    fig_removal_methods = plotComparison(data_info,data_removal_methods,'Data cleaning methods comparison');

    %% decide for removal method
    data_info = getRemovalMethod(data_info,data_removal_methods);

    %% remove non-stimuli response
    data = removeNonStimuliResponse(data,data_info.settings.non_stim_removal_method,data_info.sigWin);

    %% save data as .mat
    save([result_dir 'ProcessedData.mat'],'data');

    %% save fig_removal_methods as fig and eps and clos afterwards
    saveFigure(fig_removal_methods,[result_dir 'fig_removal_methods']);
end



function data_obj = prepareDataObject(data_info,data,ROI,number_bootstrapsamples)
    %% make data obj
    data_obj = data_handle_corrected(data_info,data,ROI);
    %% set bootstrap samples
    data_obj.prepare_samples_array(number_bootstrapsamples)
end

function [data_info,data_obj] = compareDataCleaningMethods(data_info,data_obj,result_dir,list_cleaning_methods)
    %% set default list of cleaning methods
    if nargin < 4
        list_cleaning_methods = {'none','GIF','LSM','LSM_GIF'};
    end

    %%% set result file name
    result_data_file = [ result_dir 'data_cleaning_data.mat'];

    if isfile(result_data_file)
        %% load data cleaning data
        load(result_data_file,'data_cleaning_data');
        disp('Data cleaning data already exist, loading them from file')
        disp(result_data_file)
    else
        %% calculate results for different data cleaning methods and 
        data_cleaning_data = calcDataCleaningComparisonData(data_obj, list_cleaning_methods);
        
        %% save data cleaning data
        save(result_data_file,'data_cleaning_data');
    end

    %% plot figures to compare data cleaning methods
    fig_data_cleaning = plotComparison(data_info,data_cleaning_data,'Data cleaning methods comparison');

    %% decide for cleaning method
    [data_info,data_obj] = getDataCleaningMethod(data_info,data_obj,list_cleaning_methods);

    %% save fig_data_cleaning as fig and eps and close afterwards
    saveFigure(fig_data_cleaning,[result_dir 'fig_data_cleaning']);
end

function [data_info,data_obj,spacial_analysis_results] = setFilterSettingsBasedOnSpacialAnalysis(data_info,data_obj,result_dir)
    result_file = [result_dir  'spacial_analysis_results.mat'];
    if exist(result_file,'file')
        %% load spacial anlysis results
        load(result_file,'spacial_analysis_results','data_info');
        disp('Spacial analysis results already exist, loading them from file')
        disp(result_file)

        %% plot pinwheel plateus
        fig_spacial_analysis = [];
        fig_spacial_analysis = plotSpacialAnalysis(spacial_analysis_results,fig_spacial_analysis,data_info);
    else   
        %% calculate powerspectrum/modularity/column spacing 
        [spacial_analysis_results,fig_spacial_analysis,data_info] = getSpacialAnalysis(data_info,data_obj);
        
        %% save spacial_analysis_results
        save(result_file,'spacial_analysis_results','data_info');
    end

    %% set filter settings
    [fig_spacial_analysis,data_info] = setFilterSettings(data_info,spacial_analysis_results,fig_spacial_analysis);
    data_obj.set_filter_parameters('lowpass',data_info.settings.lowpass_mm);
    data_obj.set_filter_parameters('highpass',data_info.settings.highpass_mm);

    %% save fig_spacial_analysis as fig and eps and close afterwards
    saveFigure(fig_spacial_analysis,[result_dir 'fig_spacial_analysis']);
end

function [data_info,data_obj] = refineLowpassFilterBasedOnPinwheelDensityAnalysis(data_info,data_obj,spacial_analysis_results,result_dir)
    result_file = [result_dir  'pinwheel_plateu_data.mat'];
    if exist(result_file,'file')
        %% load pinwheel plateau data
        load(result_file,'pinwheel_plateu_data','data_info');
        disp('Pinwheel plateu data already exist, loading them from file')
        disp(result_file)

        %% plot pinwheel plateus
        fig_pinwheel_plateus = [];
        fig_pinwheel_plateus = getPinwheelPlateuFigure(data_info,data_obj,pinwheel_plateu_data,spacial_analysis_results,fig_pinwheel_plateus);
    else
        %% calculate pinwheel plateus
        [pinwheel_plateu_data,fig_pinwheel_plateus,data_info] = getPinwheelPlateauData(data_info,data_obj,spacial_analysis_results);

        %% save pinwheel_plateu_data
        save([result_dir 'pinwheel_plateu_data.mat'],'pinwheel_plateu_data','data_info');
    end

    %% fine tune lowpass filter cutoff
    [data_info,fig_pinwheel_plateus] = fineTuneLowpassFilter(data_info,data_obj,pinwheel_plateu_data,spacial_analysis_results,fig_pinwheel_plateus);
    data_obj.set_filter_parameters('lowpass',data_info.settings.lowpass_mm);

    %% save fig_pinwheel_plateus as fig and eps and close afterwards
    saveFigure(fig_pinwheel_plateus,[result_dir 'fig_pinwheel_plateus']);
end

%%% get analysis data

function  data_removal_methods = calcRemovalMethodsComparisonData(data,ROI,removal_method_list,number_bootstrapsamples,alpha,filter_maps)
    %% set default parameters
    if nargin < 4
        removal_method_list = {'ff','mean','cocktail party'};
    end
    if nargin < 5
        number_bootstrapsamples = 100;
    end
    if nargin < 6
        alpha = 0.05;
    end
    if nargin < 7
        filter_maps = false;
    end

    %% calculate results for different methods
    data_removal_methods = struct();
    for ii = 1:length(removal_method_list)
        data_removal_methods.(removal_method_list{ii}) = getRemovalMethodData(data,ROI,removal_method_list{ii},number_bootstrapsamples,alpha,filter_maps);
    end

end

function removal_method_data = getRemovalMethodData(data,data_info,ROI,removal_method,number_bootstrapsamples,alpha,filter_maps)
    %% set default parameters
    if nargin < 4
        removal_method = 'cocktail party';
    end
    if nargin < 5
        number_bootstrapsamples = 100;
    end
    if nargin < 6
        alpha = 0.05;
    end
    if nargin < 7
        filter_maps = false;
    end
    %% calculate removal method data
    data = removeNonStimuliResponse(data,removal_method,data_info.sigWin);

    %% make data_obj
    data_obj = data_handle_corrected(data_info,data,ROI);

    %% set bootstrap samples
    data_obj.prepare_samples_array(number_bootstrapsamples)

    %% calculate comparison data
    removal_method_data = getComparisonData(data_obj,alpha ,filter_maps);
end

function data_cleaning_data = calcDataCleaningComparisonData(data_obj,list_cleaning_methods,alpha,filter_maps)
    %% set default parameters
    if nargin < 3
        list_cleaning_methods = {'none','GIF','LSM','LSM_GIF'};
    end
    if nargin < 4
        alpha = 0.05;
    end
    if nargin < 5
        filter_maps = false;
    end
     
    %% calculate results for different data cleaning methods
    data_cleaning_data = struct();
    for ii = 1:length(list_cleaning_methods)
        data_cleaning_data.(list_cleaning_methods{ii}) = getCleaningData(data_obj,list_cleaning_methods{ii},alpha,filter_maps);
    end
 
end
 
 
function cleaning_data = getCleaningData(data_obj,cleaning_method,alpha,filter_maps)
 
    %% set default alpha
    if nargin < 4
        alpha = 0.05;
    end
    %% set default filter maps
    if nargin < 5
        filter_maps = false;
    end

    %% set data cleaning method
    data_obj.set_data_cleaning_method(cleaning_method);

    %% calculate cleaning_data
    cleaning_data = getComparisonData(data_obj,alpha,filter_maps);
 
 
end

function comparison_data = getComparisonData(data_obj,alpha,filter_maps)

    %% save parameters
    comparison_data.alpha = alpha;
    comparison_data.filter_maps = filter_maps;

    %% calculate map
    if filter_maps
        comparison_data.z = data_obj.filter_map(data_obj.read_map());
    else
        comparison_data.z = data_obj.read_map();
    end
    comparison_data.ROI = data_obj.ROI;

    %% calculate confidence intervalls of the pereferred angle
    [comparison_data.CI_angle,~,~,BottstapSampleMaps] = getCI(data_obj,alpha,'bca',filter_maps);

    %% calculate powerprofile
    [analysis_range_mm,~,~,~]=getDefaultAnalysisRangePowerProfile(data_obj.info);
    comparison_data.power_profile = define_filter_settings(data_obj.info,comparison_data.ROI,comparison_data.z,analysis_range_mm);

    %% calculate signal to noise ratio
    comparison_data.SNR = calcSNR(BottstapSampleMaps,comparison_data.ROI,comparison_data.z);
    %comparison_data.SNR = calcSNR_OPM_Data(data_obj);
end

function SNR = calcSNR(BottstapSampleMaps,ROI,MeanMap)
    if nargin < 3
        %% calculate mean map by averaging over bootstrap samples
        MeanMap = mean(BottstapSampleMaps,3);
    end
    %% calculate DiffMaps as the difference between the bootstrap samples and the mean map
    DiffMaps = BottstapSampleMaps-MeanMap;
    %% calculate signal power as the mean of the squared mean map
    signal_power = nanmean(MeanMap(ROI).*conj(MeanMap(ROI)),'all');
    %% calculate noise power as the mean of the squared difference maps
    noise_power_2D = nanmean(DiffMaps.*conj(DiffMaps),3);
    noise_power = nanmean(noise_power_2D(ROI),'all');
    %% calculate SNR as the ratio of signal power to noise power
    SNR = signal_power/noise_power;
end

function [spacial_analysis_results,fig_spacial_analysis,data_info] = getSpacialAnalysis(data_info,data_obj)

    %% set default parameters
    [analysis_range_mm,step_size_mm,max_dim_mm,size_pix_mm]=getDefaultAnalysisRangePowerProfile(data_info);

    %% confirm or change parameters
    disp('The default analysis range for the spacial analysis is:')
    disp([num2str(min(analysis_range_mm)) ':' num2str(step_size_mm) ':' num2str(max(analysis_range_mm))])
    disp(['For comparison the largtest dimension of the map is:' num2str(max_dim_mm)])
    disp(['The size of the pixels in mm is:' num2str(size_pix_mm)])
    use_default = getBooleanInput('Do you want to use these parameters?');
    if ~use_default
        analysis_range_mm = askForAnalysisRange();
    end

    fig_spacial_analysis = [];
    confirm_results = false;

    while ~confirm_results
        spacial_analysis_results = struct();

        %% calculate power spectrum
        disp('Calculating power spectrum')
        spacial_analysis_results.power_profile = calculatePowerSpectrum(data_info,data_obj,analysis_range_mm);
        %% calculate column spacing
        disp('Calculating column spacing')
        spacial_analysis_results.column_spacing = calculateColumnSpacing(data_obj,analysis_range_mm);
        % %% calculate modularity
        % disp('Calculating modularity')
        % spacial_analysis_results.modularity = calculateModularity(data_info,data_obj,analysis_range_mm);

        %% plot results
        fig_spacial_analysis = plotSpacialAnalysis(spacial_analysis_results,fig_spacial_analysis,data_info);

        %% confirm results
        confirm_results = getBooleanInput('Do you want to use these results?');

        if ~confirm_results
            %% show previous analysis range
            disp('The previous analysis range was:')
            disp([num2str(min(analysis_range_mm)) ':' num2str(analysis_range_mm(2)-analysis_range_mm(1)) ':' num2str(max(analysis_range_mm))])
            
            %% ask for new parameters
            analysis_range_mm = askForAnalysisRange();
        end
        data_info.settings.analysis_range_mm = analysis_range_mm;
    end

    %% save settings to data_info
    data_info.settings.spacial_analysis_range_mm = analysis_range_mm;

end

function column_spacing = calculateColumnSpacing(data_obj,analysis_range_mm)
    %% convert analysis range
    smallest_w_mm = min(analysis_range_mm);
    largest_w_mm = max(analysis_range_mm);
    w_step_mm = (largest_w_mm-smallest_w_mm)/length(analysis_range_mm);

    %% get map 
    z = data_obj.read_map(); % Map

    %% calculate column spacing
    [column_spacing.average_spacing_mm,column_spacing.local_spacing_mm,column_spacing.newROI,column_spacing.WavletCoefficient] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
end

function power_spectrum = calculatePowerSpectrum(data_info,data_obj,analysis_range_mm)
    %% calculate power spectrum
    power_spectrum = define_filter_settings(data_info,data_obj.ROI,data_obj.read_map(),analysis_range_mm);
end

function [analysis_range_mm,step_size_mm,max_dim_mm,size_pix_mm]=getDefaultAnalysisRangePowerProfile(data_info)
    max_dim_mm = max(data_info.field_size_pix)/data_info.pix_per_mm;
    size_pix_mm = 1/data_info.pix_per_mm;
    step_size_mm = 0.1;
    max_range_mm = 2.5;
    if max_dim_mm < max_range_mm
        max_range_mm = max_dim_mm;
    end
    analysis_range_mm = 0.01:step_size_mm:max_range_mm;
end



function [pinwheel_plateu_data,fig_pinwheel_plateus,data_info] = getPinwheelPlateauData(data_info,data_obj,spacial_analysis_results)

    %% set default parameters
    splitROI = {};
    lowpass_cutoffs_mm = calcDefaultLowpassCutoffs(data_info,spacial_analysis_results);

    %% confirm or change default parameters
    disp('The pinwheel plateau calculation has two parameters:')
    disp('1. The lowpass cutoffs in mm')
    disp(['The default lowpass cutoffs are: ' num2str(min(lowpass_cutoffs_mm)) ':' num2str(lowpass_cutoffs_mm(2)-lowpass_cutoffs_mm(1)) ':' num2str(max(lowpass_cutoffs_mm))])
    disp('2. The original ROI can be split up in sub-ROIs, to check for local pinwheel plateaus')
    disp('This can be useful for large maps, if the pinwheel plateaus are not clearly visible.')
    disp('The default is the original ROI. If you want to split the ROI, please input the number of sub-ROIs. Or input the coordinates for each sub-ROI manually.') 
    
    confirm_parameters = getBooleanInput('Do you want to use these parameters?');
    if ~confirm_parameters
        confirm_lowpass_cutoffs = getBooleanInput('Do you want to use the default lowpass cutoffs?');
        if ~confirm_lowpass_cutoffs
            lowpass_cutoffs_mm = askLowpassRange();
        end

        confirm_splitROI = getBooleanInput('Do you want to split the ROI?');
        if ~confirm_splitROI
            splitROI = askSplitROI(data_obj);
        end
        
    end


    %% calculate pinwheel plateus
    fig_pinwheel_plateus = [];
    confirm_results = false;
    while ~confirm_results
        pinwheel_plateu_data = calculatePinwheelPlateaus(data_info,data_obj,splitROI,lowpass_cutoffs_mm);
        fig_pinwheel_plateus = getPinwheelPlateuFigure(data_info,data_obj,pinwheel_plateu_data,spacial_analysis_results,fig_pinwheel_plateus);
        confirm_results = getBooleanInput('Do you want to use these results?');

        if ~confirm_results
            %% show previous lowpass cutoffs
            disp('The previous lowpass cutoffs were:')
            disp([num2str(min(lowpass_cutoffs_mm)) ':' num2str(lowpass_cutoffs_mm(2)-lowpass_cutoffs_mm(1)) ':' num2str(max(lowpass_cutoffs_mm))])
            disp('You can change the lowpass cutoffs and/or the sub-ROIs?');
            
            %% change lowpass?
            confirm_lowpass_cutoffs = getBooleanInput('Do you want to use the previous lowpass cutoffs?');
            if ~confirm_lowpass_cutoffs
                lowpass_cutoffs_mm = askLowpassRange();
            end
            
            %% change ROI?
            confirm_splitROI = getBooleanInput('Do you want to keep the previous ROI split?');
            while ~confirm_splitROI
                splitROI = askSplitROI(data_obj);
                figure_spliROI = getSplitROIFigure(data_obj,splitROI);
                confirm_splitROI = getBooleanInput('Do you want use this split of the ROI?');
                close(figure_spliROI)
            end
        end

    end

    %% save settings to data_info
    data_info.settings.splitROI = splitROI;
    data_info.settings.lowpass_cutoffs_mm = lowpass_cutoffs_mm;

end

function DefaultLowpassCutoffs = calcDefaultLowpassCutoffs(data_info,spacial_analysis_results)
    %% set spacial scale
    spacial_scale = spacial_analysis_results.column_spacing.average_spacing_mm;

    %% set lowpass cutoffs
    DefaultLowpassCutoffs = 0.2*spacial_scale:0.01:1.1*spacial_scale;

end

function pinwheel_plateu_data = calculatePinwheelPlateaus(data_info,data_obj,splitROI,lowpass_cutoffs_mm)
    
    %% set previous filter settings
    lowpass_mm = data_info.settings.lowpass_mm;
    highpass_mm = data_info.settings.highpass_mm;

    %% calculate partial ROIs
    partialROIs = getPartialROIs(data_obj.ROI,splitROI);
    number_splitROI = size(partialROIs,3);

    %% save full ROI
    ROI = data_obj.ROI;
    
    %% calculate pinwheel plateus
    filtersPwNumbers = cell(number_splitROI,1);
    for ii = 1:number_splitROI
        data_obj.set_ROI(partialROIs(:,:,ii));
        filtersPwNumbers{ii} = find_lowpassPwNumber(data_obj,data_info,highpass_mm,lowpass_mm,lowpass_cutoffs_mm);
    end
    
    %% save results
    pinwheel_plateu_data.filtersPwNumbers = filtersPwNumbers;
    pinwheel_plateu_data.splitROI = splitROI;
    pinwheel_plateu_data.partialROIs = partialROIs;
    pinwheel_plateu_data.lowpass_cutoffs_mm = lowpass_cutoffs_mm;

    %% reset ROI in obj
    data_obj.set_ROI(ROI)
end

function partialROIs = getPartialROIs(ROI,splitROI)
    if isempty(splitROI)
        partialROIs = reshape(ROI,[size(ROI) 1]);
    else
        partialROIs = zeros(size(ROI,1),size(ROI,2),length(splitROI));
        for ii = 1:length(splitROI)
            partialROIs(:,:,ii) = calcPartialROI(splitROI{ii}.x,splitROI{ii}.y,ROI);
        end
    end
end


function splitROI = calcSplitROI(NumSplitROI,ROI,lower_limit_ROI_sizes)
    if nargin < 3
        lower_limit_ROI_sizes = 0.7;
    end
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
            if (sum(partialROI,'all')/size_splitROI) > lower_limit_ROI_sizes            
                splitROI{n}.x = splitROI_x;
                splitROI{n}.y = splitROI_y;
                n = n+1;
            end
        end
    end

end

function partROI = calcPartialROI(x,y,ROI)
    %% get Meshrid size ROI
    [X,Y] = meshgrid(1:size(ROI,2),1:size(ROI,1));
    %% get partial ROI
    partROI_bool = ((X >= x(1) & X <= x(2)) & (Y >= y(1) & Y <= y(2)));
    partROI = zeros(size(ROI));
    partROI(partROI_bool)=1;
    %% only use overlapping area with original ROI
    partROI = (partROI.*ROI);
end 

%%% input parameters
function data_info = getRemovalMethod(data_info,data_removal_methods)
    %% display input options
    disp('The options for the removal of stimulus independent signal are:')
    for ii = 1:length(data_removal_methods)
        disp([num2str(ii) '. ' data_removal_methods{ii} ': ' getDesciptionNonStimRemoval(data_removal_methods{ii})])
    end
    disp('Please choose one of the methods above.')
    %% ask for removal method
    removal_method = getNumericInput(['Input removal method (1-' num2str(length(data_removal_methods)) '): ']);
    %% check removal method
    if removal_method > length(data_removal_methods)
        error('Invalid removal method')
    end
    %% set removal method
    data_info.settings.non_stim_removal_method = data_removal_methods{removal_method};

end
function desc = getDesciptionNonStimRemoval(removal_method)
    switch removal_method
        case {'none','no','off'}
            desc = 'No removal of stimulus independent signal.';
        case {'ff','first frame','firstframe'}
            desc = 'Removes first frame for each recording.';
        case {'average','mean','average response','mean response'}
            desc = 'Removes mean over all data';
        case {'cocktail party','cocktailparty'}
            desc = 'Removes mean over all stimuli for each recording trial.';
        otherwise
            error('Invalid removal method')
    end
end

function [data_info,data_obj] = getDataCleaningMethod(data_info,data_obj,list_cleaning_methods) 

    %% describe data cleaning methods
    disp('The options for data cleaning methods are:')
    for ii = 1:length(list_cleaning_methods)
        disp([num2str(ii) '. ' list_cleaning_methods{ii} ': ' getDesciptionDataCleaning(list_cleaning_methods{ii})])
    end
    disp('Please choose one of the methods above.')
    %% ask for data cleaning method
    data_cleaning_method = getNumericInput(['Input data cleaning method (1-' num2str(length(list_cleaning_methods)) '): ']);
    %% check data cleaning method
    if data_cleaning_method > length(list_cleaning_methods)
        error('Invalid data cleaning method')
    end
    %% set data cleaning method
    data_info.settings.data_cleaning_method = list_cleaning_methods{data_cleaning_method};
    data_obj.set_data_cleaning_method(list_cleaning_methods{data_cleaning_method});
end

function desc = getDesciptionDataCleaning(cleaning_method)
    switch cleaning_method
        case 'none'
            desc = 'No data cleaning.';
        case 'GIF'
            desc = 'General Indicator Functions (GIF).';
        case 'LSM'
            desc = 'Local Similarity Minimization (LSM).';
        case 'LSM_GIF'
            desc = 'Combination of LSM and GIF.';
        otherwise
            error('Invalid data cleaning method')
    end
end


function analysis_range_mm = askForAnalysisRange()
    %% ask for analysis range
    min_analysis_range_mm = getNumericInput('Input min_analysis_range_mm: ');
    step_size_mm = getNumericInput('Input step_size_mm: ');
    max_analysis_range_mm = getNumericInput('Input max_analysis_range_mm: ');
    analysis_range_mm = min_analysis_range_mm:step_size_mm:max_analysis_range_mm;
end

function lowpass_cutoffs_mm = askLowpassRange()
    %% ask for lowpass cutoffs
    min_lowpass_cutoff_mm = getNumericInput('Input min_lowpass_cutoff_mm: ');
    step_size_mm = getNumericInput('Input step_size_mm: ');
    max_lowpass_cutoff_mm = getNumericInput('Input max_lowpass_cutoff_mm: ');
    lowpass_cutoffs_mm = min_lowpass_cutoff_mm:step_size_mm:max_lowpass_cutoff_mm;
end


function [fig_spacial_analysis,data_info] = setFilterSettings(data_info,spacial_analysis_results,fig_spacial_analysis)
    
    if doFilterSettingsExist(data_info)
        confirm_filter_settings = getBooleanInput('Should the existing filter settings be keeped?');

    else
        disp('Set filter settings based on spacial analysis results')
        disp('Lowpass filter cutoff (mm) should be smaller than the peak of the main signal')
        disp('Highpass filter cutoff (mm) should be larger than the peak of the main signal')
        confirm_filter_settings = false;
    end

    while ~confirm_filter_settings
    
        %% set filter settings
        lowpass_mm = getNumericInput('lowpass_mm: ');
        highpass_mm = getNumericInput('highpass_mm: ');
        data_info.settings.lowpass_mm = lowpass_mm;
        data_info.settings.highpass_mm = highpass_mm;

        %% plot settings
        fig_spacial_analysis = plotSpacialAnalysis(spacial_analysis_results,fig_spacial_analysis,data_info);

        %% confirm filter settings
        confirm_filter_settings = getBooleanInput('Do you want to confirm the filter settings?');
    end
end

function splitROI = askSplitROI(data_obj)
    %% display input options
    disp('You can set the split of the ROI in two ways:')
    disp('1. By inputting the number of time the ROI should be split.')
    disp('2. By inputting the coordinates for each sub-ROI manually.')
    split_ROI_input = getNumericInput('Which do you choose (1 or 2)? ');

    switch split_ROI_input
        case 1
            %% ask for number of sub-ROIs
            n_sub_ROIs = getNumericInput('Input number of ROI splits: ');
            lower_limit_ROI_sizes = 0.7;
            splitROI = calcSplitROI(n_sub_ROIs,data_obj.ROI, lower_limit_ROI_sizes);
            disp(['Final number of sub-ROIs: ' num2str(length(splitROI))])
            disp(['sub_ROIS with less than ' num2str(lower_limit_ROI_sizes*100) '% with the original ROI were removed.'])
        case 2
            %% ask for coordinates of sub-ROIs
            
        otherwise
            error('Invalid input')
    end

end

function [data_info,fig_pinwheel_plateus] = fineTuneLowpassFilter(data_info,data_obj,pinwheel_plateu_data,spacial_analysis_results,fig_pinwheel_plateus)
    disp('Fine tune lowpass filter cutoff using the pinwheel numbers')
    disp('Lowpass filter cutoff (mm) should be inside the range where pinwheel numbers are stable')
    
    if doFilterSettingsExist(data_info)
        confirm_filter_settings = getBooleanInput('Should the existing filter settings be keeped?');
    else
        confirm_filter_settings = false;
    end

    while ~confirm_filter_settings
        %% set filter settings
        lowpass_mm = getNumericInput('new lowpass cut-off [mm]: ');
        
        data_info.settings.lowpass_mm = lowpass_mm;

        %% plot setting
        fig_pinwheel_plateus = getPinwheelPlateuFigure(data_info,data_obj,pinwheel_plateu_data,spacial_analysis_results,fig_pinwheel_plateus);

        %% confirm filter settings
        confirm_filter_settings = getBooleanInput('Do you want to confirm the filter settings?');
    end

end

%%% plot results

function saveFigure(fig,filename,fig_size)
    %% set default figure size (cm)
    if nargin < 3
        fig_size = [0 0 20 40]; % [left, bottom, width, height] in cm
    end
    %% save figure as .fig 
    saveas(fig,[filename '.fig']);
    %% save figure as .eps
    fig.Units = 'centimeters';
    fig.Position = fig_size;  % [left, bottom, width, height] in cm

    fig.PaperPositionMode = 'auto';  % Keeps on-screen size when printing
    print('-depsc', '-r900', [filename '.eps']);
    %% close figure
    close(fig)
end

function fig_comparison = plotComparison(data_info,comparison_data, main_title, list_methods) 
    
    if nargin < 4
        %% extract from struct
        list_methods = fieldnames(comparison_data);
        %list_clenaing_methods =  {'none','GIF','LSM','LSM_GIF'};
    end
    
    %% prepare figure and tiled layout
    fig_comparison = figure;
    t = tiledlayout(4,length(list_methods));
    title(t,[main_title ' ' data_info.animal ' ' data_info.ID])

    %% plot data cleaning methods
    for ii = 1:length(list_methods)
        title_str = getTitleCleaningMethod(list_methods{ii});
        plotComparisonPlots(comparison_data.(list_methods{ii}),t,title_str)
    end

end

function title_str = getTitleCleaningMethod(cleaning_method)
    switch cleaning_method
        case 'none'
            title_str = 'No cleaning';
        otherwise
            title_str = replace(cleaning_method,'_',' and ');
    end
end

function plotComparisonPlots(data_cleaning_data,t,title_str)
    %% extract data from data_obj
    z = data_cleaning_data.z;
    ROI = data_cleaning_data.ROI;
    CI_angle = data_cleaning_data.CI_angle;
    alpha = data_cleaning_data.alpha;
    SNR = data_cleaning_data.SNR;
    power_profile = data_cleaning_data.power_profile;


    %% get max map range
    [x_range,y_range] = getRangeXY_ROI(ROI);

    %% plot map
    nexttile(t)
    plot_map(z,ROI,0,1)
    xlim(x_range)
    ylim(y_range)
    title(['Map ' title_str ])

    %% plot power spectrum
    nexttile(t)
    plotPowerspectrum(power_profile) 

    %% plot confidence interval
    str_alpha = num2str(100*(1-alpha));
    ax = nexttile(t);
    plot_mapAbs(CI_angle,[str_alpha '% CI angle [°]'],180,0,ROI,ax);
    xlim(x_range)
    ylim(y_range)

    %% plot SNR, mean CI and % CI < 90° as txt
    nexttile(t)
    text(0.1,0.8,['SNR: ' num2str(SNR)],'Units','normalized')
    text(0.1,0.6,['mean CI: ' num2str(mean(CI_angle(ROI)))] ,'Units','normalized')
    text(0.1,0.4,['% CI < 90°: ' num2str(sum(CI_angle(ROI) < 90)/sum(ROI(:))*100)] ,'Units','normalized')
    % remove axis 
    axis off

end



function fig_spacial_analysis = plotSpacialAnalysis(spacial_analysis_results,fig_spacial_analysis,data_info)
    %% remove old figure
    if ~isempty(fig_spacial_analysis)
        close(fig_spacial_analysis)
    end

    %% prepare figure and tiled layout
    fig_spacial_analysis = figure;
    t = tiledlayout(3,1);
    
    %% plot power spectrum
    nexttile(t)
    plotPowerspectrum(spacial_analysis_results.power_profile)
    plotFilterSettings(data_info,spacial_analysis_results.power_profile.values_kspace,true)

    nexttile(t)
    plotPowerprofile(spacial_analysis_results.power_profile) 
    plotFilterSettings(data_info,spacial_analysis_results.power_profile.values,false)

    %% plot mean wavlet 
    nexttile(t)
    plotWaveletCoeefficients(spacial_analysis_results.column_spacing)
    plotFilterSettings(data_info,spacial_analysis_results.column_spacing.WavletCoefficient.Y_mean,false)

end

function plotWaveletCoeefficients(column_spacing)
    %% get column spacing data
    average_spacing_mm = column_spacing.average_spacing_mm;
    wavelet_coefficients = column_spacing.WavletCoefficient;
    
    %% plot
    plot(wavelet_coefficients.X/1000, wavelet_coefficients.Y_mean,'DisplayName','mean wavelet');
    hold on
    plot(wavelet_coefficients.XI/1000, wavelet_coefficients.YI_mean,'DisplayName','mean wavelet fit');
    hold on
    plot([average_spacing_mm average_spacing_mm],[min(wavelet_coefficients.Y_mean) max(wavelet_coefficients.Y_mean)],':black','DisplayName','average spacing');
    xlabel('scale [mm]')
    ylabel('mean wavelet coefficient')
end

function plotPowerspectrum(power_profile) 
    plot(power_profile.k_mm_inv,power_profile.values_kspace)
    xlabel('k [1/mm]')
    ylabel('power')
    title('power spectrum')
    xlim([0 5])
end

function plotPowerprofile(power_profile) 
    plot(power_profile.scale_mm,power_profile.values)
    xlabel('scale [mm]')
    ylabel('power')
    title('power profile')
end

function plotFilterSettings(data_info,y_values,inverse)
    if nargin < 3
        inverse = false;
    end

    %% set min and max y values
    min_y = min(y_values);
    max_y = max(y_values);

    %% plot filter settings and plot them as vertical lines
    if doFilterSettingsExist(data_info)
        lowpass = data_info.settings.lowpass_mm;
        highpass = data_info.settings.highpass_mm;

        if inverse
            lowpass = 1/lowpass;
            highpass = 1/highpass;
        end

        hold on
        plot([lowpass lowpass],[min_y max_y],'-r','DisplayName','lowpass cutoff')
        hold on
        plot([highpass highpass],[min_y max_y],'-b','DisplayName','highpass cutoff')
    end
end

function existFilterSettings = doFilterSettingsExist(data_info)
    if isfield(data_info,'settings')
        if isfield(data_info.settings,'lowpass_mm') && isfield(data_info.settings,'highpass_mm')
            existFilterSettings = true;
        else
            existFilterSettings = false;
        end
    else
        existFilterSettings = false;
    end
end


function fig_pinwheel_plateus = getPinwheelPlateuFigure(data_info,data_obj,pinwheel_plateu_data,spacial_analysis_results,fig_pinwheel_plateus) 
    %% remove old figure
    if ~isempty(fig_pinwheel_plateus)
        close(fig_pinwheel_plateus)
    end

    %% prepare figure and tiled layout
    fig_pinwheel_plateus = figure;
    t = tiledlayout(3,1);
    title(t,['Refine lowpass filter cutoff ' data_info.animal ' ' data_info.ID])

    %% plot power spectrum
    nexttile(t)
    plotPowerprofile(spacial_analysis_results.power_profile)
    plotFilterSettings(data_info,spacial_analysis_results.power_profile.values,false)
    xlim([min(pinwheel_plateu_data.lowpass_cutoffs_mm) max(pinwheel_plateu_data.lowpass_cutoffs_mm)])

    %% plot pinwheel plateus
    nexttile(t)
    plotPinwheelPlateus(pinwheel_plateu_data.filtersPwNumbers)
    plotFilterSettings(data_info,pinwheel_plateu_data.filtersPwNumbers{1}.global_plateau.lowpass_vs_density(:,2),false)
    xlim([min(pinwheel_plateu_data.lowpass_cutoffs_mm) max(pinwheel_plateu_data.lowpass_cutoffs_mm)])

    %% plot split ROIs
    nexttile(t)
    plot_map(data_obj.filter_map(data_obj.read_map()),data_obj.ROI,0,1)
    plotSplitROIs(pinwheel_plateu_data.splitROI)
    [x_range,y_range] = getRangeXY_ROI(data_obj.ROI);
    xlim(x_range)
    ylim(y_range)
    
end

function figure_spliROI = getSplitROIFigure(data_obj,splitROI) 

    %% prepare figure
    figure_spliROI = figure;
    z = data_obj.filter_map(data_obj.read_map());
    [x_range,y_range] = getRangeXY_ROI(data_obj.ROI);
    ROI = data_obj.ROI;
    plot_map(z,ROI,0,1)
    hold on
    plotSplitROIs(splitROI)
    xlim(x_range)
    ylim(y_range)
end

function plotPinwheelPlateus(filtersPwNumbers)
    color_list = getColorList(length(filtersPwNumbers));
    for ii = 1:length(filtersPwNumbers)
        plot(filtersPwNumbers{ii}.global_plateau.lowpass_vs_density(:,1),filtersPwNumbers{ii}.global_plateau.lowpass_vs_density(:,2),'Color',color_list{ii})%,'DisplayName',['pinwheel data ROI ' num2str(ii)]
        hold on
    end
    title('pinwheel plateau plot')
    ylabel('pinwheel number')
    xlabel('lowpass filter cutoff [mm]')
end

function plotSplitROIs(splitROI)
    if length(splitROI) > 1
        color_list = getColorList(length(splitROI));
        for ii = 1:length(splitROI)
            rectangle('Position',[splitROI{ii}.x(1),splitROI{ii}.y(1),splitROI{ii}.x(2)-splitROI{ii}.x(1),splitROI{ii}.y(2)-splitROI{ii}.y(1)],'EdgeColor',color_list{ii},'LineWidth',2)
        end
        title('split up ROIs')
    end
    title('ROI')
end

function color_list = getColorList(n)
    if n < 7
        color_list = {'r','g','b','c','m','y'};
    else
        %% generate a list of n random colors 
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%OLD FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function data_obj = setROI(data_obj)
%     %% future features
%     % plot CI
%     % plot area that encircles areas with CI smaller x
%     %apply data cleaning?
%     %% use internal function
%     data_obj.make_ROI()
% end


% function [data_info,PwPlateus] = fineTuneLowpassFilter_old(data_info,data_obj,power_profile)
%     %% get OPM
%     z = data_obj.filter_map(data_obj.read_map());
%     ROI = data_obj.ROI;

%     %% initialize figure and parameters
%     f = figure;
%     partialROIs = [];
%     correct_PwPlateu = false;
%     lowpass_cutoffs = 0.1:0.1:2;

%     disp('For large maps it is recommended to only use one or multiple  smaller ROIs for the pinwheel plateu calculation.')
%     while ~correct_PwPlateu

%         %% if partialROIs are already set ask if they should be used
%         if ~isempty(partialROIs)
%             close(f)
%             f = plot_partialROIs(z,ROI,partialROIs);
%             redoPartialROI = getBooleanInput('Do you want to redo the partial ROIs?');
%         else
%             redoPartialROI = true;
%         end

%         if redoPartialROI
%             %% plot map and ask for partial ROI
%             close(f)
%             [partialROIs,f] = getPartialROIs(z,ROI);
%         end

%         %% check lowpass_cutoffs
%         disp('The lowpass_cutoffs to calculate the pinwheel plateus are:') 
%         disp([num2str(min(lowpass_cutoffs)) ':'  max(lowpass_cutoffs) 'the mean step size is:' num2str(mean(diff(lowpass_cutoffs)))])
%         disp(['the current lowpass_cutoff is: ' num2str(data_info.settings.lowpass_mm)])
%         correct_lowpass_cutoffs = getBooleanInput('Do you want to use these lowpass_cutoffs?');
%         if ~correct_lowpass_cutoffs
%             lowpass_cutoffs = input('Input different lowpass_cutoffs: ');
%         end

%         %% calculate pinwheel plateus
%         PwPlateus = calculatePinwheelPlateus(data_obj,partialROIs);
%         close(f)
%         f = plotPwPlateus(PwPlateus,data_info.settings.lowpass_mm);
%         correct_PwPlateu = getBooleanInput('Are the pinwheel plateus usable or do you want to redo them ?');
%     end

%     %% fine tune lowpass filter cutoff
%     disp('Fine tune lowpass filter cutoff using the pinwheel numbers')
%     disp('Lowpass filter cutoff (mm) should be inside the pinwheel plateus')
%     correctFilter = false;
%     while ~correctFilter
%         lowpass_mm = getNumericInput('lowpass_mm: ');
%         close(f)
%         f = plotPwPlateusAndPowerProfile(PwPlateus,power_profile,partialROIs,data_obj,lowpass_mm,data_info.settings.highpass_mm);
%         correctFilter = getBooleanInput('Is the lowpass cutoff correct?');
%     end
%     close(f)

%     data_info.settings.lowpass_mm = lowpass_mm;
    
% end


% function f = plotPwPlateusAndPowerProfile(PwPlateus,power_profile,partialROIs,data_obj,lowpass_mm,highpass_mm,SavFig)
%     if nargin < 7
%         SavFig = '';
%     end


%     f = figure;
%     t = tiledlayout(3,1);

%     %% plot power profile
%     nexttile(t)
%     plot(power_profile.scale_mm,power_profile.values);
%     hold on
%     plot([lowpass_mm lowpass_mm],[0 max(power_profile.values)],'r','DisplayName','lowpass cutoff')
%     hold on
%     plot([highpass_mm highpass_mm],[0 max(power_profile.values)],'r','DisplayName','highpass cutoff')
%     xlabel('Scale in mm')
%     ylabel('Power')
%     legend()
%     xlim([0 max(profile_range_mm)])
%     title('Power profile')

%     %% plot pinwheel plateus
%     color_list = getColorList(length(PwPlateus));
%     nexttile(t)
%     for ii = 1:length(PwPlateus)
%         plot(PwPlateus{ii}.global_plateau.lowpass_vs_density(:,1),PwPlateus{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)],'Color',color_list{ii})
%         hold on
%     end
%     plot([lowpass_mm lowpass_mm],[0 max(PwPlateus.values)],'r','DisplayName','current lowpass cutoff')

%     %% plot ROIs
%     nexttile(t)
%     z = data_obj.filter_map(data_obj.read_map());
%     [x_range,y_range] = getRangeXY_ROI(ROI);
%     ROI = data_obj.ROI;
%     plot_map(z,ROI,0,1)
%     for ii = 1:size(partialROIs,3)
%         contour(partialROIs(:,:,ii),[1 1],'LineWidth',2,'Color',color_list{ii});
%     end
%     xlim(x_range)
%     ylim(y_range)
%     title('ROIs')

%     if ~isempty(SavFig)
%         print(f, '-dpsc', '-fillpage', [SaveFig, '.ps']);
%         close(f);
%     end
% end

% % function color_list = getColorList(n)
% %     color_list = cell(n,1);
% %     for ii = 1:n
% %         color_list{ii} = [rand rand rand];
% %         for jj = 1:ii-1
% %             while sum(abs(color_list{ii}-color_list{jj})) < 0.5
% %                 color_list{ii} = [rand rand rand];
% %             end
% %         end
% %     end
% % end 

% function f = plotPwPlateus(PwPlateus,lowpass_cutoff)
%     f = figure;
%     for ii = 1:length(PwPlateus)
%         plot(PwPlateus{ii}.global_plateau.lowpass_vs_density(:,1),PwPlateus{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)])
%         hold on
%     end
%     plot([lowpass_cutoff lowpass_cutoff],[0 max(PwPlateus.values)],'r','DisplayName','current lowpass cutoff')
%     xlabel('lowpass cutoff in mm')
%     ylabel('density')
%     legend()
% end

% function PwPlateus = calculatePinwheelPlateus(data_obj,partialROIs,lowpass_cutoffs)

%     PwPlateus = cell(length(splitROI),1);
%     for ii = 1:length(splitROI)
%         data_obj.set_ROI(partialROIs(:,:,ii));
%         PwPlateus{ii} = find_lowpassPwNumber(data_obj,data_info,highpass_mm,lowpass_mm,lowpass_cutoffs);
%         % plot(PwPlateus{ii}.global_plateau.lowpass_vs_density(:,1),PwPlateus{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)],'Color',color_list{ii})
%         % hold on
%     end

% end


% % function [partialROIs,f] = getPartialROIs(z,ROI)
% %     f = figure;
% %     plot_map(z,ROI,0,1)
% %     doPartialROI = getBooleanInput('Do you want to use a partial ROI? ');
% %     if doPartialROI
% %         disp('Select partial ROIs by hand using the format [[x1 y1 x2 y2];[x1 y1 x2 y2]] or input a number of ROIs to select automatically!')

% %         CorrectPartialROIs = false;
% %         while ~CorrectPartialROIs
% %             %% convert partialROIs to uniform format
% %             SplitROI = getNumericInput('Input partialROIs: ');
% %             partialROIs = calculatePartialROIs(SplitROI);

% %             close(f)
% %             f = plot_partialROIs(z,ROI,partialROIs);
% %             CorrectPartialROIs = getBooleanInput('Do you want to use these partial ROIs? ');
% %         end
% %     else
% %         partialROIs = reshape(ROI,[size(ROI,1),size(ROI,2),1]);
% %     end
% % end

% function f = plot_partialROIs(z,ROI,partialROIs)
%     [x_range,y_range] = getRangeXY_ROI(ROI);
%     f = figure;
%     plot_map(z,ROI,0,1)
%     hold on
%     for ii = 1:size(partialROIs,3)
%         contour(partialROIs(:,:,ii),[1 1],'white','LineWidth',2);
%     end
%     xlim(x_range)
%     ylim(y_range)
% end

% function partialROIs = calculatePartialROIs(SplitROI)
%     if length(splitROI)>1
%         partialROIs = calcPartialROIs(data_obj.ROI,splitROI);
%     else
%         partialROIs = calcSplitROI(data_obj.ROI,SplitROI);
%     end

% end

% % function partialROIs = calcSplitROI(ROI,NumSplitROI)
% %     %% get ranges ROI
% %     [x_range,y_range] = getRangeXY_ROI(ROI);
% %     %% get split ROI
% %     x_split = round(linspace(x_range(1),x_range(2),NumSplitROI+1));
% %     y_split = round(linspace(y_range(1),y_range(2),NumSplitROI+1));
% %     size_splitROI = (x_split(2)-x_split(1))*(y_split(2)-y_split(1));
% %     n=1;
% %     for ii = 1:NumSplitROI
% %         for jj = 1:NumSplitROI
% %             %% get splitROI
% %             splitROI_x = [x_split(ii) x_split(ii+1)];
% %             splitROI_y = [y_split(jj) y_split(jj+1)];

% %             %% calc ROI
% %             partialROI = calcPartialROI(splitROI_x,splitROI_y,ROI);

% %             %% check size of splitROI
% %             if (sum(partialROI,'all')/size_splitROI) > 0.8            
% %                 % splitROI{n}.x = splitROI_x;
% %                 % splitROI{n}.y = splitROI_y;
% %                 partialROIs(:,:,n) = partialROI;
% %                 n = n+1;
% %             end
% %         end
% %     end

% % end

% % function partROI = calcPartialROI(x,y,ROI)
% %     %% get Meshrid size ROI
% %     [X,Y] = meshgrid(1:size(ROI,2),1:size(ROI,1));
% %     partROI_bool = ((X >= x(1) & X <= x(2)) & (Y >= y(1) & Y <= y(2)));
% %     partROI = zeros(size(ROI));
% %     partROI(partROI_bool)=1;
% %     partROI = (partROI.*ROI);
% % end                                                                                                                                                                                                                

% function partialROIs = calcPartialROIs(ROI,splitROI)
%     partialROIs = zeros(size(ROI,1),size(ROI,2),length(splitROI));
%     for ii = 1:size(splitROI,1)
%         partialROIs(:,:,ii) = calcPartialROI(splitROI(ii,[1 3]),splitROI(ii,[2 4]),ROI);
%     end
% end

% %contour(A, [0.5 0.5], 'r', 'LineWidth', 2);



% function [data_info,power_profile] = setFilterSettingsPowerSpectrum(data_info,data)

%     %% parameters for profile
%     max_dim_mm = max(data_info.field_size_pix)/data_info.pix_per_mm;
%     profile_range_mm = 0.1:0.01:0.8*max_dim_mm;
    
%     %% calculate power profile

%     GoodPowerprofile = false;
%     f= figure();
%     while ~GoodPowerprofile
%         close(f)
%         power_profile = define_filter_settings(data_info,ones(size(data,1:2)),data_obj.data,profile_range_mm);
%         f = plotPowerProfile(power_profile);
%         GoodPowerprofile = getBooleanInput('Do you want to use this Powerprofile?');
%         if ~GoodPowerprofile
%             profile_range_mm = input('Input different profile_range_mm (0.1:0.01:5): ');
%         end
%     end

%     disp('Set filter settings based on power profile!')
%     disp('Lowpass filter cutoff (mm) should be smaller than the peak of the main signal')
%     disp('Highpass filter cutoff (mm) should be larger than the peak of the main signal')
%     lowpass_mm = getNumericInput('lowpass_mm: ');
%     highpass_mm = getNumericInput('highpass_mm: ');
%     hold on
%     plot([lowpass_mm lowpass_mm],[0 max(power_profile.values)],'r','DisplayName','lowpass cutoff')
%     hold on
%     plot([highpass_mm highpass_mm],[0 max(power_profile.values)],'r','DisplayName','highpass cutoff')
%     legend()

%     saveFilter = false;
%     while ~saveFilter
%         saveFilter = input('Save filter settings? (y/n): ','s');
%         if saveFilter
%             data_info.settings.lowpass_mm = lowpass_mm;
%             data_info.settings.highpass_mm = highpass_mm;
%             close(f)
%         else
%             lowpass_mm = getNumericInput('lowpass_mm: ');
%             highpass_mm = getNumericInput('highpass_mm: ');
%             close(f)
%             f = plotPowerProfile(power_profile);
%             hold on
%             plot([lowpass_mm lowpass_mm],[0 max(power_profile.values)],'r','DisplayName','lowpass cutoff')
%             hold on
%             plot([highpass_mm highpass_mm],[0 max(power_profile.values)],'r','DisplayName','highpass cutoff')
%             legend()
%         end
%     end
% end

% function f = plotPowerProfile(power_profile)
%     f = figure;
%     plot(power_profile.scale_mm,power_profile.values);
%     xlabel('Scale in mm')
%     ylabel('Power')
%     xlim([0 max(profile_range_mm)])
%     set(gca,'fontsize',15)
% end



% function updateDataInfoList(data_info)
%     %% get dir of this script
%     script_dir = fileparts(mfilename('fullpath'));

%     %% pull latest version of script_dir from git
%     system(['git pull ' script_dir]);

%     %% get data_info_file_list path
%     data_info_file_list_file_path = [script_dir data_info.data_info_file_list];

%     %% read in data_info_file_list
%     data_info_file_list = readtable(data_info_file_list_file_path);

%     %% get values Dataset,MakeInfo,DestinationFolder,ExperimentID,ExperimentName
%     Dataset = data_info.animal;
%     MakeInfo = data_info.data_info_file_path;
    
%     if isfield(data_info,'DestinationFolder')
%         DestinationFolder = data_info.DestinationFolder;
%     else
%         DestinationFolder = '';
%     end

%     ExperimentID = data_info.exp_num;
%     ExperimentName = data_info.ID

%     %% check whether dataset already exists
%     if any(strcmp(data_info_file_list.Dataset,Dataset) & strcmp(data_info_file_list.ExperimentID,ExperimentID))
%         disp(['Dataset ' Dataset ' with ExperimentID ' num2str(ExperimentID) ' already exists in data_info_file_list'])
%         %% disp full MakeInfo and ExperimentName in table
%         disp(data_info_file_list(strcmp(data_info_file_list.Dataset,Dataset) & strcmp(data_info_file_list.ExperimentID,ExperimentID),:))

%     else
%         disp(['Adding new dataset ' Dataset ' with ExperimentID ' num2str(ExperimentID) ' to data_info_file_list'])
%         %% add new row to table. columns: Dataset,MakeInfo,DestinationFolder,ExperimentID,ExperimentName
%         data_info_file_list = [data_info_file_list;{Dataset,MakeInfo,DestinationFolder,ExperimentID,ExperimentName}];

%         %% write table to data_info_file_list
%         writetable(data_info_file_list,data_info_file_list_file_path);

%         %% commit changes to git
%         system(['git add ' data_info_file_list_file_path ; 'git commit -m "added new dataset enty ' Dataset num2str(ExperimentID) '"']);
%         disp('Changes saved and commited to git')
%     end  
% end