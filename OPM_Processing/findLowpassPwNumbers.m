function ResultData = findLowpassPwNumbers(data_obj,lowpass_cutoffs,density,column_spacing,sigma)
%    case {'mouse lemur'}
%        lowpass_mm = 0.2;
%         highpass_mm = 1.2;
%         lowpass_cutoffs = linspace(0.1, 0.79,70);
%% input parameter
if nargin <3
    density=false;
end

if nargin < 4
    if density == false
        column_spacing = 0;
    else
        error('To calculate pinwheel density, the columnspacing is needed!')
    end
end
if nargin < 5
    sigma = 0.1;
end

%% save lowpass parameter
lowpass_mm = data_obj.filter_parameters.lowpass;
ResultData.lowpass_cutoffs = lowpass_cutoffs;


%% get pinwheel density for filter settings
z = data_obj.read_map;
for ii = 1:length(lowpass_cutoffs)
    data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(ii))
    z_filtered =data_obj.filter_map(z);
    if density == false
        [NumberPw,~,~,~,~,~,~] = find_pinwheels(z_filtered,0,data_obj.ROI,0);
        ResultData.NumberPw(ii)=NumberPw;
    elseif density == true
        [DensityPw,NumberPw] = calcPwDensity(z_filtered,column_spacing,data_obj.ROI,data_obj.info.pix_per_mm,sigma);
        ResultData.NumberPw(ii)=NumberPw;
        ResultData.DensityPw(ii)=DensityPw;
    end
end
data_obj.set_filter_parameters('lowpass',lowpass_mm)

% if fit
%     % set axis and fit plateau to global data
%     x = lowpass_cutoffs;%.*data_info.pixels_per_mm./average_w;
%     y = pw_number;
%     
%     % figure
%     % plot(x,y)
%     
%     min_length = 0.1;%*data_info.pixels_per_mm./average_w;
%     fit_range = [min(x) max(x)];
%     
%     [pwd, plt, params, ssel, min_error]= fit_piecewise_linear_Manuel(x,y,fit_range,min_length,false);
%     
%     
%     
%     filters.global_plateau.lowpass_vs_density = [lowpass_cutoffs(:) pw_number(:)];
% end

end


function [DensityPw,NumberPw] = calcPwDensity(z,local_spacing_mm,ROI,pix_per_mm,sigma)
    [NumberPw,~,~,PWxList,PWyList,~, ~] = find_pinwheels(z,0,ROI,0);
    average_spacing_mm = mean(local_spacing_mm(ROI));
    local_pw_dens = put_gaussians(size(ROI,1),size(ROI,2), PWxList, PWyList,average_spacing_mm*pix_per_mm,sigma,ROI);
    local_pw_dens = local_pw_dens./sum(local_pw_dens(ROI)).*NumberPw;
    LocalPwDensityFixedFilter = local_pw_dens.*(local_spacing_mm*pix_per_mm).^2;
    DensityPw = mean(LocalPwDensityFixedFilter(ROI));
end
















