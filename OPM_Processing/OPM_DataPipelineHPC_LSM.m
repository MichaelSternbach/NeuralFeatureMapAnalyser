function OPM_DataPipelineHPC_LSM(animal,experiment_num,AnimalDataFolder,DataFolderMain,getCI,Bootstrapsamples,scale,setFilterParameter,ColumnSpacingSteps,PinwheelPlateuSteps_mm,RelativeSizeGaussKernel,ConfidenceLevel)
    %% inputs
    %smallest_w_mm,w_step_mm,largest_w_mm
%     experiment_num = 4;
%     animal = 'Dunnart';
%     
%     AnimalDataFolder = '~/CIDBN1/'; 
    
    %% check arg in
    disp('check arg in')

    if nargin <4
        DataFolderMain = 'Data/';
    end
    
    if nargin <5
        getCI = true;
    end
    if nargin <6
        Bootstrapsamples = 100;
    end
    if nargin <7
        scale = 0.3;
    end
    
    if nargin <8
        setFilterParameter = false;
    end

    %% parameter spacing finder
    if nargin <9
        smallest_w_mm = 0.1;
        w_step_mm = 0.05;
        largest_w_mm = 1.5;
    else
        ColumnSpacingSteps = checkFormatNum(ColumnSpacingSteps);
        smallest_w_mm = ColumnSpacingSteps(1);
        w_step_mm = ColumnSpacingSteps(2);
        largest_w_mm = ColumnSpacingSteps(3);
    end
    
    %% CI confidence parameter
    if nargin < 11
        ConfidenceLevel = 0.05; %p-value
    end
    
    %% parameter pinwheel density calculations
    
    if nargin <13
        RelativeSizeGaussKernel=0.5; %relative size gaussian kernel
    end
    
    
    %% check formats
    disp('check formats')
    experiment_num = checkFormatNum(experiment_num);
%     getCI = boolean(getCI);
    Bootstrapsamples = checkFormatNum(Bootstrapsamples);
    scale = checkFormatNum(scale);
    smallest_w_mm = checkFormatNum(smallest_w_mm);
    w_step_mm = checkFormatNum(w_step_mm);
    largest_w_mm = checkFormatNum(largest_w_mm);
    RelativeSizeGaussKernel = checkFormatNum(RelativeSizeGaussKernel);
    
    %% make dir
    disp('make dir')
    DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
    DataFolder = strrep(DataFolder,' ','_');
    disp(DataFolder)
    mkdir(DataFolder)



    %% get animal data
    disp('get animal data')
    [~,~,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder);
    data_obj.prepare_samples_array(Bootstrapsamples)
    
    
    %% deactivate LSM
    data_obj.apply_LSM(false);
    
    %% get Noise Covarienaces unfiltered, without LSM
    disp('get Noise Covarienaces unfiltered')
    DoFilter = false;
    getNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],'vector',DoFilter,scale);
    getNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],'align',DoFilter,scale);
    
    %% get Noise Covarienaces filtered, without LSM
    disp('get Noise Covarienaces unfiltered')
    DoFilter = true;
    getNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],'vector',DoFilter,scale);
    getNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],'align',DoFilter,scale);

    
    %% activate LSM
    data_obj.apply_LSM(true);
    
    %% get column spacing
    disp('get column spacing')
%     getCI = true;%true;
    [mean_spacing_mm,local_spacing_mm,newROI] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI);
    % test bootstrapping
    

    %% get CI filtered
    disp('get CI filtered')
    DoFilter = true;
    calcCIs(data_obj,ConfidenceLevel,DoFilter,DataFolder);
    
    
    %% get pinwheel infos
    if nargin <11
        PinwheelPlateuSteps_mm = linspace(0.2*mean_spacing_mm, mean_spacing_mm,50);
    else
        PinwheelPlateuSteps_mm = checkFormatNum(PinwheelPlateuSteps_mm);
    end
    disp('get pinwheel infos')
%     getCI = true;%true;
    do_plotting=0;
    PwInfo = getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,PinwheelPlateuSteps_mm,RelativeSizeGaussKernel);
    
    
%     %% get CI unfiltered
%     disp('get CI unfiltered')
%     DoFilter = false;
%     calcCIs(data_obj,alpha,DoFilter,DataFolder);
% 
    %% get Noise Covarienaces unfiltered
    disp('get Noise Covarienaces unfiltered')
    DoFilter = false;
    getNoiseCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
    getNoiseCovariances(data_obj,DataFolder,'align',DoFilter,scale);
% 
% 
    %% get Noise Covarienaces filtered
    disp('get Noise Covarienaces filtered')
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


    disp(['mean spacing [mm] ' num2str(mean_spacing_mm)])
    disp(['mean pw density ' num2str(PwInfo.MeanPwDensity)])
    disp(['mean pw number ' num2str(PwInfo.NumberPw)])


    
    disp('Finishesd!')
end































