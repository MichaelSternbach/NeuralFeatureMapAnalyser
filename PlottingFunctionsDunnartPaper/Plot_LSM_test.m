function Plot_LSM_test(animal,experiment_num,AnimalDataFolder,DataFolderMain,FigureFoler)
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

    
    DiffType = 'vector';
    
    %% check formats
    disp('check formats')
    experiment_num = checkFormatNum(experiment_num);

    
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
    [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,trial_ii,DoRectangleROI,AnimalDataFolder);
    %data_obj.prepare_samples_array(Bootstrapsamples)
    
    FigureFile = [FigureFoler 'LSM/LSM_' animal data_obj.info.ID];
    
    %% plot Noise Covarienaces unfiltered, without LSM
    disp('get Noise Covarienaces unfiltered')
    DoFilter = false;
%     getNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],'vector',DoFilter,scale);
%     getNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],'align',DoFilter,scale);
    Cov = loadNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],DiffType,DoFilter);


    

    %% get Noise Covarienaces unfiltered
    disp('get Noise Covarienaces unfiltered')
    DoFilter = false;
%     getNoiseCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
%     getNoiseCovariances(data_obj,DataFolder,'align',DoFilter,scale);
    CovLSM = loadNoiseCovariances(data_obj,[DataFolder],DiffType,DoFilter);
    
    %% plot Noise Covarienaces filtered, without LSM
    disp('get Noise Covarienaces unfiltered')
    DoFilter = true;
%     getNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],'vector',DoFilter,scale);
%     getNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],'align',DoFilter,scale);
    CovFiltered = loadNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],DiffType,DoFilter);
    
    
    %% get Noise Covarienaces filtered
    disp('get Noise Covarienaces filtered')
    DoFilter = true;
%     getNoiseCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
%     getNoiseCovariances(data_obj,DataFolder,'align',DoFilter,scale);
    CovFilteredLSM = loadNoiseCovariances(data_obj,[DataFolder 'withoutLSM'],DiffType,DoFilter);


    f = figure();
    t = tiledlayout(2,4);
    title(t,data_info.ID)

 
    
    ax=nexttile;
    plot_mapAbs(abs(Cov.Var.C1).*Cov.ROI,'No Filter & No LSM',nan,nan,data_obj.ROI,ax)
    ax=nexttile;
    plot_mapAbs(abs(CovLSM.Var.C1).*CovLSM.ROI,'No Filter & LSM applied',nan,nan,data_obj.ROI,ax)
    ax=nexttile;
    plot_mapAbs(abs(CovFiltered.Var.C1).*CovFiltered.ROI,'Filtered & No LSM',nan,nan,data_obj.ROI,ax)
    ax=nexttile;
    plot_mapAbs(abs(CovFilteredLSM.Var.C1).*CovFilteredLSM.ROI,'Filtered & LSM applied',nan,nan,data_obj.ROI,ax)
    
    nexttile;
    data_obj.apply_LSM(false);
    plot_map(data_obj.read_map(),data_obj.ROI,0,1)
    nexttile;
    data_obj.apply_LSM(true);
    plot_map(data_obj.read_map(),data_obj.ROI,0,1)
    nexttile;
    data_obj.apply_LSM(false);
    plot_map(data_obj.filter_map(data_obj.read_map()),data_obj.ROI,0,1)
    nexttile;
    data_obj.apply_LSM(true);
    plot_map(data_obj.filter_map(data_obj.read_map()),data_obj.ROI,0,1)
    
    savefig(f,[FigureFile '.fig'])
    print(f,'-depsc2', [FigureFile,'.eps']);
    
    disp('Finishesd!')
end































