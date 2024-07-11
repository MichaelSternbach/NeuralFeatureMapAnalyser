function Covariances = getNoiseCovariances(data_obj,DataFolder,DiffType,DoFilter,scale)
    if nargin < 4
        scale = .1;
    end
    
    if DoFilter
        CovariancesFile = [DataFolder DiffType 'Filtered_Covariances_' data_obj.info.ID '.mat'];
    else
        CovariancesFile = [DataFolder DiffType 'Covariances_' data_obj.info.ID '.mat'];
    end
    
    if isfile(CovariancesFile)
        load(CovariancesFile,'Covariances')
    else
        %% calc difference maps aka noise
        [DiffMaps,~,ROI]= getDifferenceMaps(data_obj,scale,DiffType,DoFilter);
        Covariances.ROI = ROI;
        Covariances.scale = convert_scale(scale,data_obj);
        
        %% calc noise variances
        Covariances.Var.C1 = mean(DiffMaps.*conj(DiffMaps),3);
        Covariances.Var.C2 = mean(DiffMaps.*DiffMaps,3); 

        %% calc distance dependent noise covariances 
        [Covariances.CoVar1D.Distances,Covariances.CoVar1D.C1,Covariances.CoVar1D.C2,Covariances.CoVar1D.NumDataList] = SpatialCovariance1D(DiffMaps,ROI);

        %% calc 2d noise covariance
        [Covariances.CoVar2D.C1,Covariances.CoVar2D.C2,Covariances.CoVar2D.N_PixelPairs] = SpatialCovariance2D(DiffMaps,ROI);
        
        %% save data
        save(CovariancesFile,'Covariances')
    end
end