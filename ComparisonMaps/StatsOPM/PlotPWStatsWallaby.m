function PlotPWStatsWallaby(BootstrapSampleFile,DataFolder,FigureFolder)
    
    addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
    addpath '/home/michael/Cloud//git/vone/MatlabCode/PatchAnalysis'


    addpath(DataFolder)%'/home/michael/Cloud/PhD/MarsupialData'
    load([DataFolder BootstrapSampleFile],'AllMaps')%'TestBig5000WallabyH.mat')

    [~,FileNameBootstrapSample,~] = fileparts(BootstrapSampleFile);
    
    IntermediatResultsFile = [FileNameBootstrapSample '_PW_Stats.mat'];
    FigureFilename = [FileNameBootstrapSample '_PW_Stats.fig'];
    
    
    if isfile([DataFolder IntermediatResultsFile])
        disp('PW Stats already exist')
        load([DataFolder IntermediatResultsFile],'pinwheel_stats','pinwheel_spurious')
    else
        tracker_obj = pinwheel_tracker;
        simple_track=true;
        [pinwheel_stats,pinwheel_spurious] = get_Wallaby_pinwheel_stats(AllMaps(:,:,1:100),tracker_obj,simple_track);
        save([DataFolder IntermediatResultsFile],'pinwheel_stats','pinwheel_spurious')
    end


    z_base = AllMaps(:,:,1);

    figure

    plot_map(z_base)

    hold
    plotPinwheelStats(pinwheel_stats,size(AllMaps,1:2))



    savefig([FigureFolder FigureFilename])
end

function plotPinwheelStats(pinwheel_stats,region)
    for i_pw = 1:getN_PW(pinwheel_stats)
        plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:),region)
    end
end

function  plotPinwheel(PWx,PWy,ProbabilityPW,region)
    ProbabilityLimitPW = .0;
    if ProbabilityPW >= ProbabilityLimitPW
        plotPosition(PWx(1),PWy(1),ProbabilityPW)
        plotConfidenceRegion(PWx,PWy,region)
    end
end

function plotPosition(PWx,PWy,ProbabilityPW)
    plot(PWx,PWy,'.white')
    text(PWx,PWy,num2str(ProbabilityPW),'Color','white')
end
function plotConfidenceRegion(PWx,PWy,region)%[156,171]
    CI = points_confidence_region(PWx,PWy,region,'hull');
    contour(CI,[1 1],'white')
end

function N_PW = getN_PW(pinwheel_stats)
    N_PW = size(pinwheel_stats.x,1);
end

