function Variances = getNoiseVariances(data_obj,DataFolder,DiffType,DoFilter,scale)
    if nargin < 5
        scale = 1;
    end
    
    if DoFilter
        CovariancesFile = [DataFolder DiffType 'Filtered_Variances_' data_obj.info.ID '.mat'];
    else
        CovariancesFile = [DataFolder DiffType 'Variances_' data_obj.info.ID '.mat'];
    end
    
    if isfile(CovariancesFile)
        load(CovariancesFile,'Variances')
    else
        %% calc difference maps
        [DiffMaps,~,ROI]= getDifferenceMaps(data_obj,scale,DiffType,DoFilter);
        Variances.ROI = ROI;
        Variances.scale = convert_scale(scale,data_obj);
        %% variances
        Variances.Var.C1 = mean(DiffMaps.*conj(DiffMaps),3);
        Variances.Var.C2 = mean(DiffMaps.*DiffMaps,3); 

%         %% calc distance dependent covariances 
%         [Variances.CoVar1D.Distances,Variances.CoVar1D.C1,Variances.CoVar1D.C2,Variances.CoVar1D.NumDataList] = SpatialCovariance1D(DiffMaps,ROI);
% 
%         %% calc 2d covariance
%         [Variances.CoVar2D.C1,Variances.CoVar2D.C2,Variances.CoVar2D.N_PixelPairs] = SpatialCovariance2D(DiffMaps,ROI);
        
        %% save data
        save(CovariancesFile,'Variances')
    end
end