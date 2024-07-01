% experiment_num = 1;%28;
% alpha=0.05;
% apply_filter = true;
% trial_ini=1;
% 
% [data_info,data_path] = info_handle('ferret',experiment_num);
% set_blocks = data_info.protocol.blocks;
% trials_to_use = find(set_blocks>0);
% 
% 
% data_obj = data_handle_corrected(data_info,[data_path,'Processed_2/trial_',num2str(trials_to_use(trial_ini)),'.mat'],[data_path,'exp_info.mat']);
% 
% stats_ini = load([data_path,'Analyzed_2/characterization/',sprintf('trial_%d_domain_stats',trials_to_use(trial_ini)),'.mat'],'orientation_stats','parameters');
% data_obj.set_samples_array(stats_ini.parameters.samples_array);
% data_obj.apply_LSM
% 
% orientation_stats = get_orientation_stats(data_obj,alpha,apply_filter);
%save('orientationStatsFerret.mat','orientation_stats')

load('orientationStatsFerret.mat');


%%% Plot results
addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'

%%% Plot results
addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
addpath ~/Cloud/git/vone/MatlabCode/PatchAnalysis

figure();
tiledlayout(3,2)
%tiledlayout(4,2)

nexttile
plot_map(orientation_stats(:,:,2))
title('orientation preference map')

nexttile
Abs = abs(orientation_stats(:,:,2))./mean(abs(orientation_stats(:,:,2)),'all');
plot_mapAbs(Abs,'selectivity [<selectivity>]',max(Abs,[],'all'),min(Abs,[],'all'))


nexttile
CI_Abs = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,1)))./mean(abs(orientation_stats(:,:,2)),'all');
plot_mapAbs(CI_Abs,'uncertainty selectivity [<selectivity>]',max(CI_Abs,[],'all'),min(CI_Abs,[],'all'))

% % nexttile
% % CI_Abs = abs(abs(orientation_stats(:,:,1))-abs(orientation_stats(:,:,2)))./mean(abs(orientation_stats(:,:,2)),'all');
% % plot_mapAbs(CI_Abs,'uncertainty selectivity (upper CI) [<selectivity>]',max(CI_Abs,[],'all'),min(CI_Abs,[],'all'))
% % 
% % nexttile
% % CI_Abs = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,2)))./mean(abs(orientation_stats(:,:,2)),'all');
% % plot_mapAbs(CI_Abs,'uncertainty selectivity (lower CI) [<selectivity>]',max(CI_Abs,[],'all'),min(CI_Abs,[],'all'))

% % figure();
% % tiledlayout(2,2)

% nexttile
% CI_Abs_low = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,2)))./mean(abs(orientation_stats(:,:,2)));
% plot_mapAbs(CI_Abs_low,'uncertainty selectivity (lower CI) [<selectivity>]',max(CI_Abs_low,[],'all'),min(CI_Abs_low,[],'all'))


nexttile
CI = abs(angle(orientation_stats(:,:,3)./orientation_stats(:,:,1)))/pi*90;
plot_mapAbs(CI,'uncertainty orientation [°]',90,0)

% % nexttile
% % CI = abs(angle(orientation_stats(:,:,1)./orientation_stats(:,:,2)))/pi*90;
% % plot_mapAbs(CI,'uncertainty orientation (upper CI) [°]',90,0)
% % 
% % nexttile
% % CI = abs(angle(orientation_stats(:,:,3)./orientation_stats(:,:,2)))/pi*90;
% % plot_mapAbs(CI,'uncertainty orientation (lower CI) [°]',90,0)


nexttile
plotContourAngleDelta(10,orientation_stats(:,:,2),CI)

nexttile
plotContourAngleDelta(20,orientation_stats(:,:,2),CI)