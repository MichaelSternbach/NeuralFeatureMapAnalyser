function Covariances = loadNoiseCovariances(data_obj,DataFolder,DiffType,DoFilter)

    if DoFilter
        CovariancesFile = [DataFolder DiffType 'Filtered_Covariances_' data_obj.info.ID '.mat'];
    else
        CovariancesFile = [DataFolder DiffType 'Covariances_' data_obj.info.ID '.mat'];
    end
    
    load(CovariancesFile,'Covariances')

end