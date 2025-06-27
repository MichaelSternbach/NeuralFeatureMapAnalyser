function AnalysisPipelineOPM(data_obj,DataFolder,Confidence,getCI_PwCS,Bootstrapsamples,scale)

    %% default arguments
    if nargin < 3
        Confidence = 0.05;
    end
    if nargin < 4
        getCI_PwCS = true;
    end
    if nargin < 5
        Bootstrapsamples = 400;
    end
    if nargin < 6
        scale = 100;
    end

    %% set bootstrapsamples numbers
    if isstruct(Bootstrapsamples)
        BS_PwCs = Bootstrapsamples.BS_PwCs;
        BS_HypTest = Bootstrapsamples.BS_HypTest;
        BS_CoVar = Bootstrapsamples.BS_CoVar;
    else
        BS_PwCs = Bootstrapsamples;
        BS_HypTest = Bootstrapsamples;
        BS_CoVar = Bootstrapsamples;
    end


    %% display settings

    %% analyse OPM properties and their confidence
    [PwInfo,~,~,~] = analysePinwheelsColumnSpacingCI(data_obj,DataFolder,Confidence,getCI_PwCS,BS_PwCs,SizeGaussKernelPwDensityCalc);


    %% test modularity and pinwheel arrangement
    if ~isnan(BS_HypTest)
        runHypothesisTests(data_obj,DataFolder,Bootstrapsamples,PwInfo);
    end
    
    %% analyse noise covariance
    if ~isnan(BS_CoVar)
        analyseNoiseCovariance(data_obj,DataFolder,scale,BS_CoVar);
    end

end

function [PwInfo,mean_spacing_mm,local_spacing_mm,newROI] = analysePinwheelsColumnSpacingCI(data_obj,DataFolder,Confidence,getCI_PwCS,Bootstrapsamples,SizeGaussKernelPwDensityCalc)
    
    %% set bootstrapsamples Pinwdensity
    data_obj.prepare_samples_array(Bootstrapsamples);
    disp(['BS ' num2str(size(data_obj.samples_array,3))])

    %% set analysis range
    spacial_analysis_range_mm = data_info.settings.spacial_analysis_range_mm;
    smallest_w_mm = min(spacial_analysis_range_mm);
    largest_w_mm = max(spacial_analysis_range_mm);
    w_step_mm = (largest_w_mm-smallest_w_mm)/length(spacial_analysis_range_mm);


    %% get column spacing filtered
    disp('get column spacing filtered')
    [mean_spacing_mm,local_spacing_mm,newROI] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI_PwCS,true);
    
    
    %% set lowpass cutoffs
    lowpass_cutoffs_mm= data_info.settings.lowpass_cutoffs_mm;

    %% get pinwheel infos
    disp('get pinwheel infos')
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
    do_plotting=0;
    PwInfo = getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI_PwCS,do_plotting,lowpass_cutoffs_mm,SizeGaussKernelPwDensityCalc);

    %% get CI filtered
    disp('get CI filtered')
    data_obj.prepare_samples_array(BS_CI);
    disp(['BS ' num2str(size(data_obj.samples_array,3))])
    DoFilter = true;
    calcCIs(data_obj,Confidence,DoFilter,DataFolder);

    %% get CI unfiltered
    disp('get CI unfiltered')
    DoFilter = false;
    calcCIs(data_obj,Confidence,DoFilter,DataFolder);


    %% print some results
    disp(['mean spacing [mm] ' num2str(mean_spacing_mm)])
    disp(['mean pw density ' num2str(PwInfo.MeanPwDensity)])
    disp(['mean pw number ' num2str(PwInfo.NumberPw)])
    
    %% make plots
    plotPwCsResults(data_info,data_obj,DataFolder,getCI)
end


function runHypothesisTests(data_obj,DataFolder,Bootstrapsamples,PwInfo)

    %% set analysis range
    analysis_range_mm = data_info.settings.spacial_analysis_range_mm;

    %% set bootstrapsamples for Modularity and Pinwheel tests
    data_obj.prepare_samples_array(Bootstrapsamples);
    disp(['BS ' num2str(size(data_obj.samples_array,3))])

    %% calculate column spacing unfiltered
    disp('calculate column spacing unfiltered')
    getCI = false;
    smallest_w_mm = min(analysis_range_mm);
    largest_w_mm = max(analysis_range_mm);
    w_step_mm = (largest_w_mm-smallest_w_mm)/length(analysis_range_mm);
    [mean_spacing_mm,~,~] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,true);

    %% testModularityOPM
    disp('testModularityOPM')
    testModularityOPM(data_obj,DataFolder,mean_spacing_mm,analysis_range_mm,Bootstrapsamples)


    %% testPWsOPM
    disp('testPWsOPM')
    testPWsOPM(data_obj,PwInfo.pinwheel_stats,Bootstrapsamples,DataFolder)

    %% make plots
    % ...
end


function analyseNoiseCovariance(data_obj,DataFolder,scale,Bootstrapsamples)
        %% prepare bootstrapsamples for Covariances
        data_obj.prepare_samples_array(Bootstrapsamples); 
        disp(['BS ' num2str(size(data_obj.samples_array,3))])

        %% get Noise Covarienaces unfiltered
        disp('get Noise Covarienaces unfiltered')
        DoFilter = false;
        getNoiseCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
        getNoiseCovariances(data_obj,DataFolder,'align',DoFilter,scale);
  
        %% get Noise Covarienaces filtered
        disp('get Noise Covarienaces filtered')
        DoFilter = true;
        getNoiseCovariances(data_obj,DataFolder,'vector',DoFilter,scale);
        getNoiseCovariances(data_obj,DataFolder,'align',DoFilter,scale);

        %% make plots
        % ...
end

function plotPwCsResults(data_info,data_obj,DataFolder,getCI)
        FigureFile = [DataFolder 'Factsheet_' data_info.animal ' ' data_info.ID];
        disp(FigureFile)
        if getCI
            PlotFactSheetPage(data_info,data_obj,DataFolder,FigureFile)
        else
            PlotResultWoCI(animal,data_info,data_obj,DataFolder,FigureFile)
        end
end