function OPM_DataPipelineHPC_Faster(animal,experiment_num,AnimalDataFolder,DataFolderMain,getCI,getCor,Bootstrapsamples,DataCleaning,scale,setFilterParameter,...
    ColumnSpacingCalcSteps,PwDensitCalcSteps,Confidence,SizeGaussKernelPwDensityCalc)
%     OPM_DataPipelineHPC_Faster('cat','1','~/CIDBN/','~/Test/','0','100','none','100','false',...
%     '0.1|0.05|1.5','0.2|0.02|1.2','0.05','0.5')
    %% check arg 
    disp('check arg in')

    % experiment number
    experiment_num = checkFormatNum(experiment_num);

    % main data folder
    if nargin <4
        DataFolderMain = 'Data/';
    end

    % defines whether to get confidence intervals (boolean)
    if nargin <5
        getCI = true;
    else
        getCI = checkFormatNum(getCI);
        getCI = (getCI ==1);
    end

    % defines whether to get correlations (boolean)
    if nargin <6
        getCor = true;
    else
        getCor = checkFormatNum(getCor);
        getCor = (getCor ==1);
    end
    
    % number of bootstrapsamples
    if nargin <7
        Bootstrapsamples = 100;
    else
        Bootstrapsamples = convertChar(Bootstrapsamples);
    end

    % Data cleaning method
    if nargin <8
        DataCleaning = 'none';
    end

    % determines size of the maps for the noise covariance calculation
    % if scale<1, the size is determined by the scale factor times the original map size
    % if scale>1, the size is in pixels
    if nargin <9
        scale = 100;
    else
        scale = checkFormatNum(scale);
    end

    % set filter parameter
    if nargin <10
        setFilterParameter = false;
    else
        setFilterParameter = convertChar(setFilterParameter);
        if ~isstruct(setFilterParameter)
            setFilterParameter = false;
        end
    end
    
    %% parameter spacing finder - input format: [start, step, end] in mm
    if nargin <11
        smallest_w_mm = 0.1;
        w_step_mm = 0.05;
        largest_w_mm = 1.5;
    else
        ColumnSpacingCalcSteps = checkFormatNumList(ColumnSpacingCalcSteps);
        smallest_w_mm = ColumnSpacingCalcSteps(1);
        w_step_mm = ColumnSpacingCalcSteps(2);
        largest_w_mm = ColumnSpacingCalcSteps(3);
    end

    %% length scales parameter pinwheel density calculations - input format: [start, step, end] in mm
    if nargin <12
        PwDensitCalcSteps = linspace(0.2, 1.2,50);
    else
        PwDensitCalcSteps = checkFormatNumList(PwDensitCalcSteps);
        PwDensitCalcSteps = PwDensitCalcSteps(1):PwDensitCalcSteps(2):PwDensitCalcSteps(3);
    end
    %% CI confidence parameter (p-value)
    if nargin < 13
        Confidence = 0.05;
    else
        Confidence = checkFormatNum(Confidence);
    end
    
    %% parameter pinwheel density calculations
    
    if nargin < 14
        SizeGaussKernelPwDensityCalc=0.5;
    else
        SizeGaussKernelPwDensityCalc = checkFormatNum(SizeGaussKernelPwDensityCalc);
    end
    
    %% disp iput
    disp(['animal: ' animal])
    disp(['experiment_num: ' num2str(experiment_num)])
    disp(['AnimalDataFolder: ' AnimalDataFolder])
    disp(['DataFolderMain: ' DataFolderMain])
    disp(['getCI: ' num2str(getCI)])
    disp(['getCor: ' num2str(getCor)])
    disp(['Bootstrapsamples: ' jsonencode(Bootstrapsamples)])
    disp(['DataCleaning: ' DataCleaning])
    disp(['scale: ' num2str(scale)])
    disp(['setFilterParameter: ' jsonencode(setFilterParameter)])
    disp(['ColumnSpacingCalcSteps: ' num2str(smallest_w_mm) ':' num2str(w_step_mm) ':' num2str(largest_w_mm)])
    disp(['PwDensitCalcSteps: ' num2str(PwDensitCalcSteps(1)) ':' num2str(PwDensitCalcSteps(2)-PwDensitCalcSteps(1)) ':' num2str(PwDensitCalcSteps(end))])
    disp(['SizeGaussKernelPwDensityCalc: ' num2str(SizeGaussKernelPwDensityCalc)])
    disp(['Confidence: ' num2str(Confidence)])

    


    %% set BS numbers
    if isstruct(Bootstrapsamples)
        BS_ModTest = Bootstrapsamples.BS_ModTest;   
        BS_PwTest = Bootstrapsamples.BS_PwTest;
        BS_PwDens = Bootstrapsamples.BS_PwDens;
        BS_CI = Bootstrapsamples.BS_CI;
        if getCor
            BS_Cov = Bootstrapsamples.BS_Cov;
        end
    else
        BS_ModTest = Bootstrapsamples;
        BS_PwTest = Bootstrapsamples;
        BS_PwDens = Bootstrapsamples;
        BS_Cov = Bootstrapsamples;
        BS_CI = Bootstrapsamples;
    end

    
    
    %% make dir
    disp('make dir')
    DataFolder = DataFolderMain;%[DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
    DataFolder = strrep(DataFolder,' ','_');
    disp(DataFolder)
    mkdir(DataFolder)

    %% set rng
    seed = 1234;
    rng(seed)

    %% get animal data
    disp('get animal data')
    [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder);

    %% set noise reduction algorithm
    switch lower(DataCleaning)
        case 'none'
            data_obj.apply_LSM(false);
        case 'lsm'
            data_obj.apply_LSM(true);
        case 'gif'
            data_obj.generateCleanedDataSamplesGIF()
        case 'gif_jk'
            data_obj.generateCleanedDataSamplesGIF_JK()
        otherwise
            disp(['Value DataCleaning: ' DataCleaning'])
            error('DataCleaning not recognized')
    end
    

    %% set bootstrapsamples Pinwdensity
    data_obj.prepare_samples_array(BS_PwDens);


    %% set filter parameter
    if setFilterParameter && ~isstruct(setFilterParameter)
        data_obj.set_filter_parameters('lowpass',setFilterParameter.lowpass_mm)
        data_obj.set_filter_parameters('highpass',setFilterParameter.highpass_mm)
    end

    %% get column spacing
    disp('get column spacing')
    disp(getCI)
    if getCI
        disp(['BS ' num2str(size(data_obj.samples_array,3))])
    end
    [mean_spacing_mm,local_spacing_mm,newROI] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI);
    % test bootstrapping
    
    
    %% get pinwheel infos
    disp('get pinwheel infos')
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
    do_plotting=0;
    PwInfo = getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,PwDensitCalcSteps,SizeGaussKernelPwDensityCalc);
    
    
    %% get CI filtered
    disp('get CI filtered')
    data_obj.prepare_samples_array(BS_CI);
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
    DoFilter = true;
    calcCIs(data_obj,Confidence,DoFilter,DataFolder);

    %     %% get CI unfiltered
%     disp('get CI unfiltered')
%     DoFilter = false;
%     calcCIs(data_obj,alpha,DoFilter,DataFolder);


    %% testModularityOPM
    disp('testModularityOPM')
    data_obj.prepare_samples_array(BS_ModTest);
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
    profile_range_mm = smallest_w_mm:w_step_mm:largest_w_mm;
    testModularityOPM(data_obj,DataFolder,mean_spacing_mm,profile_range_mm,BS_ModTest)
    
    %% testPWsOPM
    disp('testPWsOPM')
    data_obj.prepare_samples_array(BS_PwTest);
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
    testPWsOPM(data_obj,PwInfo.pinwheel_stats,PwInfo.pinwheel_spurious,BS_PwTest,DataFolder)

    



    if getCor
        %% prepare bootstrapsamples for Covariances
        data_obj.prepare_samples_array(BS_Cov); 

        %% get Noise Covarienaces unfiltered
        disp('get Noise Covarienaces unfiltered')
        disp(['BS ' num2str(size(data_obj.samples_array,3))])
        DoFilter = false;
        getNoiseCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
        getNoiseCovariances(data_obj,DataFolder,'align',DoFilter,scale);
  
        %% get Noise Covarienaces filtered
        disp('get Noise Covarienaces filtered')
        disp(['BS ' num2str(size(data_obj.samples_array,3))])
        DoFilter = true;
        getNoiseCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
        getNoiseCovariances(data_obj,DataFolder,'align',DoFilter,scale);

    %     %% get Map Covarienaces unfiltered
    %     disp('get Map Covarienaces unfiltered')
    %     DoFilter = false;
    %     getMapCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
    %     getMapCovariances(data_obj,DataFolder,'align',DoFilter,scale);
    % 
    %     %% get Map Covarienaces filtered
    %     disp('get Map Covarienaces filtered')
    %     DoFilter = true;
    %     getMapCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
    %     getMapCovariances(data_obj,DataFolder,'align',DoFilter,scale);
    end


    disp(['mean spacing [mm] ' num2str(mean_spacing_mm)])
    disp(['mean pw density ' num2str(PwInfo.MeanPwDensity)])
    disp(['mean pw number ' num2str(PwInfo.NumberPw)])


    
    disp('Finished!')

    disp('plot results')
    PlotPwCI(animal,data_info,data_obj,DataFolder)
end































