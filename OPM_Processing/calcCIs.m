function CI = calcCIs(data_obj,alpha,DoFilter,DataFolder)  
    if DoFilter
        CIFile = [DataFolder 'Filtered_CI_' data_obj.info.ID '.mat'];
    else
        CIFile = [DataFolder 'CI_' data_obj.info.ID '.mat'];
    end
    
    if isfile(CIFile)
        load(CIFile,'CI')
    else
        %% calc CI maps
        [CI.BCA.CI_angle,CI.BCA.CI_Abs,CI.BCA.ROI] = getCI(data_obj,alpha,'bca',DoFilter);
        [CI.SE.CI_angle,CI.SE.CI_Abs,CI.SE.ROI] = getCI(data_obj,alpha,'se',DoFilter);
        %% save data
        save(CIFile,'CI')
    end
end