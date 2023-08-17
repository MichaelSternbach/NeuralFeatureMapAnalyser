function [DiffMaps,MeanMap,ROI]= getDifferenceMaps(data_obj,scale,DiffType)   
    if nargin == 1
        scale = 1;
        DiffType = 'vector';
    elseif nargin == 2
        DiffType = 'vector';
    end
    
    [BottstapSampleMaps,MeanMap,ROI]= getBootstrapSampleMaps(data_obj,scale);
    
    switch lower(DiffType)
        case{'angle','orientation'}
            DiffMaps = angle(BottstapSampleMaps./MeanMap);
        case{'abs','selectivity'}
            DiffMaps = abs(BottstapSampleMaps)-abs(MeanMap);
        case{'vector','complex'}
            DiffMaps = BottstapSampleMaps-MeanMap;
        case{'rotated','aligned'}
            DiffMaps = (BottstapSampleMaps-MeanMap)./(MeanMap./abs(MeanMap));
    end
end
