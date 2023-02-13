% addpath '/home/michael/Cloud/PhD/MarsupialData'
% load('TestBig5000WallabyH.mat')
% load('JackknifeSamplesWallabyH.mat')
% 
% 
% orientation_stats = get_orientation_statsWallabyJason(AllMaps,JackknifeData);
% save('orientationStatsWallaby.mat','orientation_stats')

load('orientationStatsWallaby.mat');

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