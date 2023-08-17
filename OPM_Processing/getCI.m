function [CI_angle,CI_Abs,ROI] = getCI(data_obj,alpha,method)
    if nargin == 1
        alpha = 0.05;
    end
    scale = 1;
    switch lower(method)
        case{'bca','bootstrap biase corrected','biase corrected'}
            apply_filter = true;
            orientation_stats = get_orientation_stats(data_obj,alpha,apply_filter);

            CI_Abs = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,1)))/mean(abs(orientation_stats(:,:,2)),'all');            
            CI_angle = abs(angle(orientation_stats(:,:,3)./orientation_stats(:,:,2))-angle(orientation_stats(:,:,1)./orientation_stats(:,:,2)))/pi*90;
            ROI = data_obj.ROI;
        case{'se','bootstrap standard error'}
            Z = getZ(alpha);
            
            [DiffMaps,MeanMap,ROI]= getDifferenceMaps(data_obj,scale,'abs');
            CI_Abs = (mean(DiffMaps.^2,3)).^0.5*Z*2/mean(abs(MeanMap),'all');
            
            [DiffMaps,~,~]= getDifferenceMaps(data_obj,scale,'angle');
            CI_angle = (mean(DiffMaps.^2,3)).^0.5*Z*2;
            CI_angle(find(CI_angle>2*pi))=2*pi;
            CI_angle=CI_angle/pi*90;
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
    end
end