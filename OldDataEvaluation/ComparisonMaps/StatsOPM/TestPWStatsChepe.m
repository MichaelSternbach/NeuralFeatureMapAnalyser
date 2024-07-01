function TestPWStatsChepe(animal,experiment_num,data,data_info,data_path,Bootstrapsamples,DataFolder,FigureFolder)
    addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
    addpath '/home/michael/Cloud//git/vone/MatlabCode/PatchAnalysis'
    addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/AGWolfOPMDataPipeline/'

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

    IntermediatResultsFile = [animal num2str(experiment_num) '_PW_Stats.mat'];
    FigureFilename = [animal num2str(experiment_num) '_PW_Stats.fig'];
    
    
    if isfile([DataFolder IntermediatResultsFile])
        disp('PW Stats already exist')
        load([DataFolder IntermediatResultsFile],'pinwheel_stats','data_obj')
    else
        tracker_obj = pinwheel_tracker;
        simple_track=true;
        [pinwheel_stats,pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
        save([DataFolder IntermediatResultsFile],'pinwheel_stats','pinwheel_spurious','data_obj')
    end
    
    
    base = 1;
    z_base = data_obj.filter_map(data_obj.read_map(base));

    figure

    plot_map(z_base.*data_obj.ROI)

    hold
    plotPinwheelStats(pinwheel_stats,data_info.field_size_pix)
    
    savefig([FigureFolder FigureFilename])

end
    


% stats_ini = load([data_path,'Analyzed_2/characterization/',sprintf('trial_%d_domain_stats',trials_to_use(trial_ini)),'.mat'],'orientation_stats','parameters'); 
% experiment_num = 1;%28;
% animal = 'ferret';
% alpha=0.05;
% apply_filter = true;
% trial_ini=1;
% simple_track = true;
% 
% [data_info,data_path] = info_handle(animal,experiment_num);
% 
% set_blocks = data_info.protocol.blocks;
% trials_to_use = find(set_blocks>0);
% data = load([data_path,'Processed_2/trial_',num2str(trials_to_use(trial_ini)),'.mat'],'data');

function plotPinwheelStats(pinwheel_stats,field_size_pix)
    for i_pw = 1:getN_PW(pinwheel_stats)
        plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:),field_size_pix)
    end
end

function  plotPinwheel(PWx,PWy,ProbabilityPW,field_size_pix)
    ProbabilityLimitPW = .50;
    if ProbabilityPW >= ProbabilityLimitPW
        plotPosition(PWx(1),PWy(1),ProbabilityPW)
        plotConfidenceRegion(PWx,PWy,field_size_pix)
    end
end

function plotPosition(PWx,PWy,ProbabilityPW)
    plot(PWx,PWy,'.white')
    text(PWx,PWy,num2str(ProbabilityPW),'Color','white')
end
function plotConfidenceRegion(PWx,PWy,field_size_pix)
    CI = points_confidence_region(PWx,PWy,field_size_pix,'hull');
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
