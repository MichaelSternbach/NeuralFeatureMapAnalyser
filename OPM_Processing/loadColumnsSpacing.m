function [average_spacing_mm,local_spacing_mm,newROI,CI_average_spacing_mm,CI_local_spacing_mm,WavletCoefficient] = loadColumnsSpacing(data_obj,DataFolder,getCI,getWavletCoefficient)
    

    %% get mean spacing
    SpacingFile = [DataFolder 'MapSpacing_' data_obj.info.ID '.mat'];
    if getWavletCoefficient==false
        load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI')
    else
        load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI','WavletCoefficient')
    end
    if nargin < 3
        getCI = false;
    end
    %% get CIs
   if getCI == true
        CISpacingFile = [DataFolder 'CI_MapSpacing_' data_obj.info.ID '.mat'];
        load(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm','average_spacings_mm','local_spacings_mm','newROIs')
    end
end