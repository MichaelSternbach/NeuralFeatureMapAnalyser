function Covariances = getMapCovariances(data_obj,DataFolder,DiffType,DoFilter,scale)
    if nargin < 4
        scale = .1;
    end
    
    if DoFilter
        CovariancesFile = [DataFolder DiffType 'Filtered_CovariancesMaps_' data_obj.info.ID '.mat'];
    else
        CovariancesFile = [DataFolder DiffType 'CovariancesMaps_' data_obj.info.ID '.mat'];
    end
    
    if isfile(CovariancesFile)
        load(CovariancesFile,'Covariances')
    else
        %% calc bootstrap sample maps
        [Maps,~,ROI]= getMaps(data_obj,scale,DiffType,DoFilter);
        Covariances.ROI = ROI;
        Covariances.scale = convert_scale(scale,data_obj);
        
        %% variances
        Covariances.Var.C1 = mean(Maps.*conj(Maps),3);
        Covariances.Var.C2 = mean(Maps.*Maps,3); 

        %% calc distance dependent covariances 
        [Covariances.CoVar1D.Distances,Covariances.CoVar1D.C1,Covariances.CoVar1D.C2,Covariances.CoVar1D.NumDataList] = SpatialCovariance1D(Maps,ROI);

        %% calc 2d covariance
        [Covariances.CoVar2D.C1,Covariances.CoVar2D.C2,Covariances.CoVar2D.N_PixelPairs] = SpatialCovariance2D(Maps,ROI);
        
        %% save data
        save(CovariancesFile,'Covariances')
    end
end