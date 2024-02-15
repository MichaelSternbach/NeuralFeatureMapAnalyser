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
    
    
    %% check formats
    disp('check formats')
    experiment_num = checkFormatNum(experiment_num);
%     getCI = boolean(getCI);
    Bootstrapsamples = checkFormatNum(Bootstrapsamples);
    scale = checkFormatNum(scale);
    smallest_w_mm = checkFormatNum(smallest_w_mm);
    w_step_mm = checkFormatNum(w_step_mm);
    largest_w_mm = checkFormatNum(largest_w_mm);
    beta = checkFormatNum(beta);
    
    %% make dir
    disp('make dir')
    DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
    DataFolder = strrep(DataFolder,' ','_');
    disp(DataFolder)
    mkdir(DataFolder)



    %% get animal data
    disp('get animal data')
    trial_ii = 1;
    DoRectangleROI = false;
    [~,~,data_obj,~,~] = getAnimalData(animal,experiment_num,trial_ii,DoRectangleROI,AnimalDataFolder);
    data_obj.prepare_samples_array(Bootstrapsamples)

    %% get column spacing
    disp('get column spacing')
%     getCI = true;%true;
    [mean_spacing_mm,local_spacing_mm,newROI] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI);
    % test bootstrapping
    

    %% get CI filtered
    disp('get CI filtered')
    DoFilter = true;
    calcCIs(data_obj,alpha,DoFilter,DataFolder);
    
    
    %% get pinwheel infos
    if nargin <11
        llp_cutoffs = linspace(0.2*mean_spacing_mm, mean_spacing_mm,50);
    else
        llp_cutoffs = checkFormatNum(llp_cutoffs);
    end
    disp('get pinwheel infos')
%     getCI = true;%true;
    do_plotting=0;
    PwInfo = getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta);
    
    
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































