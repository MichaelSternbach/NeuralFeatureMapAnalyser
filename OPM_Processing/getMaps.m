function [Maps,MeanMap,ROI]= getMaps(data_obj,scale,DiffType,DoFilter)   
    if nargin == 1
        scale = 1;
        DiffType = 'vector';
    elseif nargin == 2
        DiffType = 'vector';
        DoFilter = true;
    elseif nargin == 3
        DoFilter = true;
    end
    
    [BottstapSampleMaps,MeanMap,ROI]= getBootstrapSampleMaps(data_obj,scale,DoFilter);
    
    switch lower(DiffType)
        case{'angle','orientation'}
            Maps = angle(BottstapSampleMaps);
        case{'abs','selectivity'}
            Maps = abs(BottstapSampleMaps);
        case{'vector','complex'}
            Maps = BottstapSampleMaps;
        case{'rotated','aligned','align'}
            Maps = (BottstapSampleMaps)./(MeanMap./abs(MeanMap));
    end
end
