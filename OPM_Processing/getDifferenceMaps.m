function [DiffMaps,MeanMap,ROI,BottstapSampleMaps]= getDifferenceMaps(data_obj,scale,DiffType,DoFilter,direction_map)
% get deviations of the bootstrap samples from the mean map
% can be used to approximate the noise and variability in the recording
    if nargin == 1
        scale = 1;
        DiffType = 'vector';
    elseif nargin == 2
        DiffType = 'vector';
        DoFilter = true;
    elseif nargin == 3
        DoFilter = true;
    end

    if nargin < 5
        direction_map = false;
    end
    
    [BottstapSampleMaps,MeanMap,ROI]= getBootstrapSampleMaps(data_obj,scale,DoFilter,direction_map);
    
    switch lower(DiffType)
        case{'angle','orientation'}
            DiffMaps = angle(BottstapSampleMaps./MeanMap);
        case{'abs','selectivity'}
            DiffMaps = abs(BottstapSampleMaps)-abs(MeanMap);
        case{'vector','complex'}
            DiffMaps = BottstapSampleMaps-MeanMap;
        case{'rotated','aligned','align'}
            DiffMaps = (BottstapSampleMaps-MeanMap)./(MeanMap./abs(MeanMap));
    end
end
