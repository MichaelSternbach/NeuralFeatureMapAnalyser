function get_filter_settings(job_ID)

if isdeployed
    
    % display host name
    disp(['Job running on ',getenv('HOSTNAME')])
    
    % save cache in scratch.local if possible
    if exist('/scratch.local/chepe/','dir')
        setenv('MCR_CACHE_ROOT','/scratch.local/chepe/')
    else
        setenv('MCR_CACHE_ROOT','/scratch04/chepe/')
        disp('No /chepe folder in scratch.local, saving in scratch04.')
    end
    
    % check that the /pairing folder can be reached
    if ~exist('/pairing/token','file')
        error('Can not access data from this node!')
    end
    
    % get job ID
    if nargin==0
        job_ID = getenv('SGE_TASK_ID');
    end
    if ischar(job_ID)
        job_ID = str2double(job_ID);
    end
end

[data_info,data_path] = zoo_keeper(job_ID);
filters = [];

%% get lowpass from design analysis

disp('Obtaining standard lowpass from design analysis...')

switch lower(data_info.animal)
    case {'aotus'}        
        lowpass_mm = 0.28;
        highpass_mm = 1.14;
        lowpass_cutoffs = linspace(0.1,0.8,71);
    case {'cat','cats'}
        lowpass_mm = 0.3;
        highpass_mm = 1.6;
        lowpass_cutoffs = linspace(0.3,1.19,90);
    case {'ferret','ferrets'}
        lowpass_mm = 0.4;
        highpass_mm = 1.5;
        lowpass_cutoffs = linspace(0.1,0.99,90);
    case {'galago','galagos'}
        lowpass_mm = 0.2;
        highpass_mm = 1.2;
        lowpass_cutoffs = linspace(0.1, 0.79,70);
    case {'shrews','shrew','treeshrews','treeshrew'}
        lowpass_mm = 0.2;
        highpass_mm = 1.2;
        lowpass_cutoffs = linspace(0.1, 0.79,70);
    case {'macaque','macaques'}
        lowpass_mm = 0.28;
        highpass_mm = 1.14;      
        lowpass_cutoffs = linspace(0.1,0.8,71);
   case {'mouse lemur'}
       lowpass_mm = 0.2;
        highpass_mm = 1.2;
        lowpass_cutoffs = linspace(0.1, 0.79,70);
end

filters.design_analysis.lowpass = lowpass_mm;
filters.design_analysis.highpass = highpass_mm;
filters.design_analysis.pw_density = data_info.design.pw_dens;

%% get lowpass with estimated plateau in global density

disp('Obtaining lowpass with global plateau fitting...')

% get pinwheel density for filter settings
pw_density_global = zeros(size(lowpass_cutoffs));
for ii = 1:length(lowpass_cutoffs)
    z_filtered = filter_map(data_info.design.z_raw,data_info.design.roi_raw==1,data_info.design.measure,lowpass_cutoffs(ii),highpass_mm);
    pw_density_global(ii) = get_pinwheel_density(z_filtered,data_info.design.local_w);
end

% set axis and fit plateau to global data
x = lowpass_cutoffs.*data_info.design.measure./data_info.design.average_w;
y = pw_density_global;

min_length = 0.15*data_info.design.measure./data_info.design.average_w;
fit_range = [0.2 1.0];

[pw_dens, center_x, plateau] = fit_piecewise_linear(x,y,min_length,fit_range);

global_lp = center_x/(data_info.design.measure./data_info.design.average_w);
global_plt = plateau/(data_info.design.measure./data_info.design.average_w);

% add to data
filters.global_plateau.lowpass = global_lp;
filters.global_plateau.highpass = highpass_mm;

filters.global_plateau.lowpass_vs_density = [lowpass_cutoffs(:) pw_density_global(:)];

filters.global_plateau.pw_density = pw_dens;
filters.global_plateau.plateau = global_plt;

%{
 % plot global plateau and fit
 clf
 plot(filters.global_plateau.plateau,filters.global_plateau.pw_density*[1 1])
 hold on
 plot(filters.global_plateau.lowpass_vs_density(:,1),filters.global_plateau.lowpass_vs_density(:,2))
 xlabel('lowpass mm')
 ylabel('pinwheel density')
%}

%% get lowpass value that gives a fixed pinwheel density (pi)

disp('Obtaining lowpass that gives pi density...')

% interpolate value that gets pi pinwheel density
cutoffs_extended = linspace(min(lowpass_cutoffs),max(lowpass_cutoffs),1000);
density_extended = interp1(lowpass_cutoffs,pw_density_global,cutoffs_extended);

[~,idx] = min(abs(density_extended-pi));
pi_lp = cutoffs_extended(idx);

% add to data
filters.pi_density.lowpass = pi_lp;
filters.pi_density.highpass = highpass_mm;
filters.pi_density.pw_density = pi; % by definition

%% get lowpass value with most pixels in corresponding plateau range

disp('Obtaining lowpass with most pixels inside plateau...')

% get range for each pixel
tmp1 = data_info.design.cutoff_lowpass(:,:,1);
tmp2 = data_info.design.cutoff_lowpass(:,:,2);
pixel_plateau = [tmp1(data_info.ROI) tmp2(data_info.ROI)];

% count fraction and find maximum
lowpass_list = linspace(0,1,100);
fraction_in_plateau = zeros(size(lowpass_list));
for ii=1:length(lowpass_list)
    fraction_in_plateau(ii) = mean(pixel_plateau(:,1)<= lowpass_list(ii) & pixel_plateau(:,2)>= lowpass_list(ii));
end
plateau_lp = median(lowpass_list(fraction_in_plateau==max(fraction_in_plateau)));

[~,idx] = min(abs(cutoffs_extended - plateau_lp));
plateau_density = density_extended(idx);

% add to data
filters.pixel_plateau.lowpass = plateau_lp;
filters.pixel_plateau.highpass = highpass_mm;
filters.pixel_plateau.pixel_plateau = pixel_plateau;
filters.pixel_plateau.pw_density = plateau_density;
filters.pixel_plateau.lowpass_vs_fraction = [lowpass_list(:) fraction_in_plateau(:)];

%{
 % plot fraction of pixels in plateau
 clf
 plot(filters.pixel_plateau.lowpass_vs_fraction(:,1),filters.pixel_plateau.lowpass_vs_fraction(:,2))
 xlabel('lowpass mm')
 ylabel('fraction in plateau')
%}

%% save results

disp('Saving results...')

pixels_per_hc = data_info.design.local_w;
save([data_path,data_info.ID,'.mat'],'filters','pixels_per_hc','-append')

end