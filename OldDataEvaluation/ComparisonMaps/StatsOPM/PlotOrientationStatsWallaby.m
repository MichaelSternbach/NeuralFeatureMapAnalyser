function PlotOrientationStatsWallaby(BootstrapSampleFile,JackknifeSampleFile,DataFolder,FigureFolder)
    addpath(DataFolder)%'/home/michael/Cloud/PhD/MarsupialData'
    load([DataFolder BootstrapSampleFile],'AllMaps')%'TestBig5000WallabyH.mat')
    load([DataFolder JackknifeSampleFile],'JackknifeData')%'JackknifeSamplesWallabyH.mat')


    [~,FileNameBootstrapSample,~] = fileparts(BootstrapSampleFile);
    
    IntermediatResultsFile = [FileNameBootstrapSample '_OrientationStats.mat'];
    FigureFilename = [FileNameBootstrapSample '_OrientationStats.fig'];
    
    
    
    if isfile([DataFolder IntermediatResultsFile])
        disp('Orientation Stats already exist')
        load([DataFolder IntermediatResultsFile],'orientation_stats')
    else
        orientation_stats = get_orientation_statsWallabyJason(AllMaps,JackknifeData);
        save([DataFolder IntermediatResultsFile],'orientation_stats')
    end
  

    %%% Plot results
    addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
    addpath ~/Cloud/git/vone/MatlabCode/PatchAnalysis

    figure();
    tiledlayout(3,2)

    nexttile
    plot_map(orientation_stats(:,:,2))
    title('orientation preference map')

    nexttile
    Abs = abs(orientation_stats(:,:,2))./mean(abs(orientation_stats(:,:,2)),'all');
    plot_mapAbs(Abs,'selectivity [<selectivity>]',max(Abs,[],'all'),min(Abs,[],'all'))


    nexttile
    CI_Abs = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,1)))./mean(abs(orientation_stats(:,:,2)),'all');
    plot_mapAbs(CI_Abs,'uncertainty selectivity [<selectivity>]',max(CI_Abs,[],'all'),min(CI_Abs,[],'all'))

    nexttile
    CI = abs(angle(orientation_stats(:,:,3)./orientation_stats(:,:,1)))/pi*90;
    plot_mapAbs(CI,'uncertainty orientation [°]',90,0)


    nexttile
    plotContourAngleDelta(10,orientation_stats(:,:,2),CI)

    nexttile
    plotContourAngleDelta(20,orientation_stats(:,:,2),CI)

    savefig([FigureFolder FigureFilename])
end