function TestOrientationStatsChepe(animal,experiment_num,data,data_info,data_path,Bootstrapsamples,DataFolder,FigureFolder)
    addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
    addpath '/home/michael/Cloud//git/vone/MatlabCode/PatchAnalysis'
    addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/AGWolfOPMDataPipeline/'
    
    alpha=0.05;
    apply_filter = true;
    
    data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
    if isa(Bootstrapsamples,'int')||(isempty(Bootstrapsamples))||isa(Bootstrapsamples,'double')
        if (isempty(Bootstrapsamples))
            Bootstrapsamples = 100;
        end
        data_obj.prepare_samples_array(Bootstrapsamples)
    else
        data_obj.set_samples_array(Bootstrapsamples.parameters.samples_array);
    end
       
    data_obj.apply_LSM

    IntermediatResultsFile = [animal num2str(experiment_num) '_OrientationStats.mat'];
    FigureFilename = [animal num2str(experiment_num) '_OrientationStats.fig'];
    
    
    if isfile([DataFolder IntermediatResultsFile])
        disp('PW Stats already exist')
        load([DataFolder IntermediatResultsFile],'orientation_stats')
    else
        orientation_stats = get_orientation_stats(data_obj,alpha,apply_filter);
        save([DataFolder IntermediatResultsFile],'orientation_stats','data_obj')
    end

    %%% Plot results

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

