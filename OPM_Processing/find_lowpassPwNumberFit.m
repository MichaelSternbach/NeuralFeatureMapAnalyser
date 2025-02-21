function [filters,plt, params, ssel] = find_lowpassPwNumberFit(data_obj,highpass_mm,lowpass_mm,lowpass_cutoffs,density,column_spacing,sigma)
%    case {'mouse lemur'}
%        lowpass_mm = 0.2;
%         highpass_mm = 1.2;
%         lowpass_cutoffs = linspace(0.1, 0.79,70);

if nargin <6
    density=false;
end

if nargin < 7
    if density == false
        column_spacing = 0;
    else
        error('To calculate pinwheel density, the columnspacing is needed!')
    end
end

if nargin < 8
    sigma = 0.1;
end

filters.design_analysis.lowpass = lowpass_mm;
filters.design_analysis.highpass = highpass_mm;

%% get lowpass with estimated plateau in global density
disp('Obtaining lowpass with global plateau fitting...')

%% get pinwheel density for filter settings
pw_number = zeros(size(lowpass_cutoffs));
z = data_obj.read_map;
for ii = 1:length(lowpass_cutoffs)
    data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(ii))
    z_filtered =data_obj.filter_map(z);
    if density == false
        [count,~,~,~,~,~,~] = find_pinwheels(z_filtered,0,data_obj.ROI,0);
        pw_number(ii)=count;
    elseif density == true
        pw_density = calcPwDensity(z_filtered,column_spacing,data_obj.ROI,data_obj.info.pix_per_mm,sigma);
        pw_number(ii)=pw_density;
    end
end
data_obj.set_filter_parameters('lowpass',lowpass_mm)
% set axis and fit plateau to global data
x = lowpass_cutoffs;%.*data_info.pixels_per_mm./average_w;
y = pw_number;

% figure
% plot(x,y)

min_length = 0.1;%*data_info.pixels_per_mm./average_w;
fit_range = [min(x) max(x)];

%[pw_dens, center_x, plateau] 
[pwd, plt, params, ssel, min_error]= fit_piecewise_linear_Manuel(x,y,fit_range,min_length,false);

% global_lp = center_x;%/(data_info.pixels_per_mm./average_w);
% global_plt = plateau;%/(data_info.pixels_per_mm./average_w);

% % add to data
% filters.global_plateau.lowpass = global_lp;
% filters.global_plateau.highpass = highpass_mm;

filters.global_plateau.lowpass_vs_density = [lowpass_cutoffs(:) pw_number(:)];

% filters.global_plateau.pw_density = pw_dens;
% filters.global_plateau.plateau = global_plt;
end


function pw_density = calcPwDensity(z,local_spacing_mm,ROI,pix_per_mm,sigma)
    [NumberPw,~,~,PWxList,PWyList,~, ~] = find_pinwheels(z,0,ROI,0);
    average_spacing_mm = mean(local_spacing_mm(ROI));
    local_pw_dens = put_gaussians(size(ROI,1),size(ROI,2), PWxList, PWyList,average_spacing_mm*pix_per_mm,sigma,ROI);
    local_pw_dens = local_pw_dens./sum(local_pw_dens(ROI)).*NumberPw;
    LocalPwDensityFixedFilter = local_pw_dens.*(local_spacing_mm*pix_per_mm).^2;
    pw_density = mean(LocalPwDensityFixedFilter(ROI));
end
















