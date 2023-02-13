addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
addpath '/home/michael/Cloud//git/vone/MatlabCode/PatchAnalysis'
addpath '/home/michael/Cloud/PhD/MarsupialData'

load('TestBig5000WallabyH.mat')


% 
% tracker_obj = pinwheel_tracker;
% simple_track=true;
%  
% [pinwheel_stats,pinwheel_spurious] = get_Wallaby_pinwheel_stats(AllMaps(:,:,1:100),tracker_obj,simple_track);
%save('WallabyStatsPW.mat','pinwheel_stats','pinwheel_spurious')

load('WallabyStatsPW.mat')

z_base = AllMaps(:,:,1);

figure

plot_map(z_base)

hold
plotPinwheelStats(pinwheel_stats)

figure
plot(sort(pinwheel_stats.probability),(1:getN_PW(pinwheel_stats))./getN_PW(pinwheel_stats))
%cumpdf(pinwheel_stats.probability)



function plotPinwheelStats(pinwheel_stats)
    for i_pw = 1:getN_PW(pinwheel_stats)
        plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:))
    end
end

function  plotPinwheel(PWx,PWy,ProbabilityPW)
    ProbabilityLimitPW = .0;
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
    CI = points_confidence_region(PWx,PWy,[156,171],'hull');
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