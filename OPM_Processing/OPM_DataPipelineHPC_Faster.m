function OPM_DataPipelineHPC_Faster(animal,experiment_num,AnimalDataFolder,DataFolderMain,getCI,Bootstrapsamples,scale,smallest_w_mm,w_step_mm,...
    largest_w_mm,llp_cutoffs,beta)
    %% inputs
    
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
    
    %% parameter spacing finder
    if nargin <10
        smallest_w_mm = 0.1;
        w_step_mm = 0.05;
        largest_w_mm = 1.5;
    end
    
    %% CI confidence parameter
    if nargin < 11
        alpha = 0.05;
    end
    
    %% parameter pinwheel density calculations
    
    if nargin <13
        beta=0.5;
    end
    
    %% disp iput
    disp(animal)
    disp(experiment_num)
    disp(AnimalDataFolder)
    disp(DataFolderMain)
    disp(getCI)
    disp(Bootstrapsamples)
    
    %% check formats
    disp('check formats')
    experiment_num = checkFormatNum(experiment_num);
    getCI = checkFormatNum(getCI);
    %getCI = boolean(getCI);
    Bootstrapsamples = convertChar(Bootstrapsamples);
    scale = checkFormatNum(scale);
    smallest_w_mm = checkFormatNum(smallest_w_mm);
    w_step_mm = checkFormatNum(w_step_mm);
    largest_w_mm = checkFormatNum(largest_w_mm);
    beta = checkFormatNum(beta);


    %% set BS numbers
    if isstruct(Bootstrapsamples)
        BS_ModTest = Bootstrapsamples.BS_ModTest;   
        BS_PwTest = Bootstrapsamples.BS_PwTest;
        BS_PwDens = Bootstrapsamples.BS_PwDens;
        BS_Cov = Bootstrapsamples.BS_Cov;
        BS_CI = Bootstrapsamples.BS_CI;
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
    [~,~,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder);

    %% set bootstrapsamples Pinwdensity
    data_obj.prepare_samples_array(BS_PwDens);

    %% get column spacing
    disp('get column spacing')
    if getCI
        disp(['BS ' num2str(size(data_obj.samples_array,3))])
    end
    [mean_spacing_mm,local_spacing_mm,newROI] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI);
    % test bootstrapping
    
    
    %% get pinwheel infos
    if nargin <11
        llp_cutoffs = linspace(0.2*mean_spacing_mm, mean_spacing_mm,50);
    else
        llp_cutoffs = checkFormatNum(llp_cutoffs);
    end
    disp('get pinwheel infos')
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
%     getCI = true;%true;
    do_plotting=0;
    PwInfo = getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta);
    
    
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

    %% get CI filtered
    disp('get CI filtered')
    data_obj.prepare_samples_array(BS_CI);
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
    DoFilter = true;
    calcCIs(data_obj,alpha,DoFilter,DataFolder);


%     %% get CI unfiltered
%     disp('get CI unfiltered')
%     DoFilter = false;
%     calcCIs(data_obj,alpha,DoFilter,DataFolder);
% 

    %% prepare bootstrapsamples for Covariances
    data_obj.prepare_samples_array(BS_Cov); 

    %% get Noise Covarienaces unfiltered
    disp('get Noise Covarienaces unfiltered')
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
    DoFilter = false;
    getNoiseCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
    getNoiseCovariances(data_obj,DataFolder,'align',DoFilter,scale);
% 
% 
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


    disp(['mean spacing [mm] ' num2str(mean_spacing_mm)])
    disp(['mean pw density ' num2str(PwInfo.MeanPwDensity)])
    disp(['mean pw number ' num2str(PwInfo.NumberPw)])


    
    disp('Finishesd!')
end































