function [CI_angle,CI_Abs,ROI,BottstapSampleMaps] = getCI(data_obj,alpha,method,apply_filter,direction_map,parallelize)
    if nargin == 1
        alpha = 0.05;
    end
    scale = 1;
    if nargin < 4
        apply_filter = true;
    end

    if nargin < 5
        direction_map = false;
    end
    if nargin < 6
        parallelize = false;
    end
    
    switch lower(method)
        case{'bca','bootstrap biase corrected','biase corrected'}
            
            [orientation_stats,BottstapSampleMaps] = get_orientation_stats(data_obj,alpha,apply_filter,direction_map,parallelize);

            CI_Abs = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,1)));%/mean(abs(orientation_stats(:,:,2)),'all')            
            
             CI_angle = abs(angle(orientation_stats(:,:,3)./orientation_stats(:,:,2))-angle(orientation_stats(:,:,1)./orientation_stats(:,:,2)));
            

            ROI = data_obj.ROI;
        case{'se','bootstrap standard error'}
            Z = getZ(alpha);
            
            %% get Diff maps
            [BottstapSampleMaps,MeanMap,ROI]= getBootstrapSampleMaps(data_obj,scale,apply_filter,direction_map,parallelize);
            DiffMapsAngle = angle(BottstapSampleMaps./MeanMap);
            DiffMapsAbs = abs(BottstapSampleMaps)-abs(MeanMap);
                    

            %[DiffMaps,MeanMap,ROI,BottstapSampleMaps]= getDifferenceMaps(data_obj,scale,'abs',apply_filter,direction_map,parallelize);
            CI_Abs = (mean(DiffMapsAbs.^2,3)).^0.5*Z*2;%/mean(abs(MeanMap),'all')
            
            %[DiffMaps,~,~]= getDifferenceMaps(data_obj,scale,'angle',apply_filter,parallelize);
            CI_angle = (nanmean(DiffMapsAngle.^2,3)).^0.5*Z*2;
            CI_angle(find(CI_angle>2*pi))=2*pi;

    end

    if direction_map
       CI_angle=CI_angle./pi*180;
    else
         CI_angle=CI_angle./pi*90;
    end

end

function Z = getZ(alpha) %https://www.mathsisfun.com/data/confidence-interval.html
    switch 1-alpha
        case .95
            Z = 1.960;
        case .99
            Z = 2.576;  
        otherwise
            disp(1-alpha)
            error('value not defined!')
    end
end