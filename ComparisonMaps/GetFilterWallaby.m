% addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
% addpath '/home/michael/Cloud//git/vone/MatlabCode/PatchAnalysis'
% addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data'/ComparisonMaps/AGWolfOPMDataPipeline/
% 
% data_set = 'Wallaby';
% experiment_num = 1;%28;
% alpha=0.05;
% apply_filter = true;
% trial_ini=1;
% simple_track = true;
% base = 1;
% 
% [data_info,data_path] = info_handle(data_set,experiment_num);
% load([data_path,data_info.ID,'.mat'],'data');
% 
% ROI = true(size(data,1:2));
% save([data_path,'/exp_info.mat'],'ROI','-append')
% 
% stim_order = data_info.stim_order;
% power_profile = define_filter_settings(data_info,data_path,data,stim_order);
map = make_map(data(:,:,:,1:40),stim_order,ROI,true);

figure
plot_map(map.*ROI)

% data_obj = data_handle_corrected(data_info,[data_path,data_info.ID,'.mat'],[data_path,'exp_info.mat']);
% 
% %stats_ini = load([data_path,'Analyzed_2/characterization/',sprintf('trial_%d_domain_stats',trials_to_use(trial_ini)),'.mat'],'orientation_stats','parameters');
% % data_obj.prepare_samples_array(100)%set_samples_array(stats_ini.parameters.samples_array);
% 
% data_obj.apply_LSM
% 
% ROI=ones(data_info.field_size_pix);
% 
% define_filter_settingsWallaby(data_set,experiment_num,data_obj)

% highpass_mm = 0.48;
% lowpass_mm = 0.1;
% 
% smallest_w = 0.2;
% largest_w = 1; 
% w_step = 0.2;
% 
% [average_w_re, local_w_re, ~] = wavelength_estimator_data(real(filter_map(data_obj.read_map,data_info.pixels_per_mm,lowpass_mm,highpass_mm)), ROI,smallest_w*1000,largest_w*1000,w_step*1000,1000/data_info.pixels_per_mm);
% disp([average_w_re local_w_re])

% 
% average_w = 171/2; 
% local_w = ROI * average_w;
% 
% lowpass_cutoffs = 0.02:0.01:0.2;
% 
% filters = find_lowpass(data_obj,data_info,highpass_mm,lowpass_mm,lowpass_cutoffs,average_w,local_w);
% 
% z_base = data_obj.filter_map(data_obj.read_map(base));
% 
% figure
% 
% plot_map(z_base)

% 
% tracker_obj = pinwheel_tracker;
% 
% [pinwheel_stats,pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
% %save('FerretStatsPW.mat','pinwheel_stats','pinwheel_spurious','data_obj')
% 
% load('FerretStatsPW.mat')
% 
% base = 1;
% z_base = data_obj.filter_map(data_obj.read_map(base));
% 
% figure
% 
% plot_map(z_base)
% 
% hold
% plotPinwheelStats(pinwheel_stats)
% 
% figure
% plot(sort(pinwheel_stats.probability),(1:getN_PW(pinwheel_stats))./getN_PW(pinwheel_stats))
% %cumpdf(pinwheel_stats.probability)



function plotPinwheelStats(pinwheel_stats)
    for i_pw = 1:getN_PW(pinwheel_stats)
        plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:))
    end
end

function  plotPinwheel(PWx,PWy,ProbabilityPW)
    ProbabilityLimitPW = .50;
    if ProbabilityPW >= ProbabilityLimitPW
        plotPosition(PWx(1),PWy(1),ProbabilityPW)
        plotConfidenceRegion(PWx,PWy)
    end
end

function plotPosition(PWx,PWy,ProbabilityPW)
    plot(PWx,PWy,'.white')
    text(PWx,PWy,num2str(ProbabilityPW),'Color','white')
end
function plotConfidenceRegion(PWx,PWy)
    CI = points_confidence_region(PWx,PWy,510,'hull');
    contour(CI,[1 1],'white')
    %plotCI(CI)
end

function N_PW = getN_PW(pinwheel_stats)
    N_PW = size(pinwheel_stats.x,1);
end

function plotCI(CI)
    for i_x = 1:size(CI,1)
        for i_y = 1:size(CI,2)
            if CI(i_x,i_y)==1
                plot(i_x,i_y,'.white')
            end
        end
    end
            
     
end
% % orientation_stats = get_orientation_stats(data_obj,alpha,apply_filter);
% % save('orientationStatsFerret.mat','orientation_stats')

% load('orientationStatsFerret.mat');

% 
% %%% Plot results
% addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
% 
% figure();
% tiledlayout(2,2)
% 
% nexttile
% plot_map(orientation_stats(:,:,2))
% title('orientation preference map')
% 
% nexttile
% Abs = abs(orientation_stats(:,:,2))./mean(abs(orientation_stats(:,:,2)),'all');
% plot_mapAbs(Abs,'selectivity [<selectivity>]',max(Abs,[],'all'),min(Abs,[],'all'))
% 
% 
% nexttile
% CI_Abs_up = abs(abs(orientation_stats(:,:,1))-abs(orientation_stats(:,:,2)))./mean(abs(orientation_stats(:,:,2)),'all');
% plot_mapAbs(CI_Abs_up,'uncertainty selectivity (upper CI) [<selectivity>]',max(CI_Abs_up,[],'all'),min(CI_Abs_up,[],'all'))
% 
% nexttile
% CI_Abs_up = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,2)))./mean(abs(orientation_stats(:,:,2)),'all');
% plot_mapAbs(CI_Abs_up,'uncertainty selectivity (lower CI) [<selectivity>]',max(CI_Abs_up,[],'all'),min(CI_Abs_up,[],'all'))
% 
% figure();
% tiledlayout(2,2)
% 
% % nexttile
% % CI_Abs_low = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,2)))./mean(abs(orientation_stats(:,:,2)));
% % plot_mapAbs(CI_Abs_low,'uncertainty selectivity (lower CI) [<selectivity>]',max(CI_Abs_low,[],'all'),min(CI_Abs_low,[],'all'))
% 
% nexttile
% CI_up = abs(angle(orientation_stats(:,:,1)./orientation_stats(:,:,2)))/pi*90;
% plot_mapAbs(CI_up,'uncertainty orientation (upper CI) [°]',90,0)
% 
% nexttile
% CI_low = abs(angle(orientation_stats(:,:,3)./orientation_stats(:,:,2)))/pi*90;
% plot_mapAbs(CI_low,'uncertainty orientation (lower CI) [°]',90,0)
% 
% 
% nexttile
% plotContourAngleDelta(10,orientation_stats(:,:,2),CI_up)
% 
% nexttile
% plotContourAngleDelta(20,orientation_stats(:,:,2),CI_up)