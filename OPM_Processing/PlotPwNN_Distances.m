function PlotPwNN_Distances(animal,experiment_Num,AnimalDataFolder,DataFolderMain,FigureFolder)
    close all
    
    %% parameter
    experiment_Num = checkFormatNum(experiment_Num);
    
%     PwDensityType= 'MeanPwDensity';
    
    %% figure folder
    mkdir(FigureFolder)
    
    FigureFile1 = [FigureFolder animal '_DistancesNN'];
    
    
% 
%     %% make Figure
% 
%     f1 = figure();
%     ax1= axes('Parent',f1);
%     
%     f2 = figure();
%     ax2= axes('Parent',f2);
%     
%     f3 = figure();
%     ax3= axes('Parent',f3);
    
    
    %% Distance Lists
    d = [];
    d_eq = [];
    d_op= [];
    all_areas = [];
    all_n = [];
    
    for experiment_num = 1: experiment_Num
        
        %% data folder
        DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
        
        %% animal
        [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,1,false,AnimalDataFolder);
        [average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,DataFolder,false);
%         z = data_obj.filter_map(data_obj.read_map());
        data_info.ID = replace(data_info.ID,'_',' ');

        %% load pinwheel data

        getCI = false;
        do_plotting=0;
        llp_cutoffs = linspace(0.01, 1,100);
        beta=0.5;

        Bootstrapsamples = 100;
        data_obj.prepare_samples_array(Bootstrapsamples)
        PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta);
        
        scale = average_spacing_mm*data_info.pix_per_mm;
        
        d = [d PwInfo.d./scale];
        d_eq = [d_eq PwInfo.d_eq./scale];
        d_op = [d_op PwInfo.d_op./scale];
        
        all_areas = [all_areas PwInfo.circ_areas]; %./scale^2
        all_n = [all_n PwInfo.n];
%         %% load CI spacing data
%         CISpacingFile = [DataFolder 'CI_MapSpacing_' data_obj.info.ID '.mat'];
%         load(CISpacingFile,'CI_average_spacing_mm')
%         
% 
%         %% load pinwheel CI data
%         CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
%         load(CIPwFile,'alpha','PwInfosBS','PwInfosJS')


        
    end
    
    plot_nn_dist('d', d, 'dunnart_nn_dist_d')
    plot_nn_dist('d++', d_eq, 'dunnart_nn_dist_deq')
    plot_nn_dist('d+-', d_op, 'dunnart_nn_dist_dop')
    plot_count_sd(all_areas, all_n,'dunnart_count_sd')
    
%     %% mouse lemur paper data
%     load microcebus_Huber.mat
% 
%     plot_nn_dist('d', d, 'microcebus_nn_dist_d')
%     plot_nn_dist('d++', d_eq, 'microcebus_nn_dist_deq')
%     plot_nn_dist('d+-', d_op, 'microcebus_nn_dist_dop')
%     plot_count_sd(all_areas, all_n,'microcebus_count_sd')
% 
%     load macaque_Angelucci.mat
% 
%     plot_nn_dist('d', d, 'macaque_nn_dist_d')
%     plot_nn_dist('d++', d_eq, 'macaque_nn_dist_deq')
%     plot_nn_dist('d+-', d_op, 'macaque_nn_dist_dop')
%     plot_count_sd(all_areas, all_n,'macaque_count_sd')
% 
%     tmp = load('macaque_GrinvaldOkamoto.mat');
%     d = [d; tmp.d];
%     d_eq = [d_eq; tmp.d_eq];
%     d_op = [d_op; tmp.d_op];
%     all_areas = [all_areas; tmp.all_areas];
%     all_n = [all_n; tmp.all_n];
% 
%     plot_nn_dist('d', d, 'macaqueAll_nn_dist_d')
%     plot_nn_dist('d++', d_eq, 'macaqueAll_nn_dist_deq')
%     plot_nn_dist('d+-', d_op, 'macaqueAll_nn_dist_dop')
%     plot_count_sd(all_areas, all_n,'macaqueAll_count_sd')
    
%     %% save figure1
%     figure(f1)
%     histogram(ax1,d,bins,'DisplayStyle','stairs','Normalization','pdf','DisplayName',data_info.ID)
%     
%     
%     xlabel('distance [Lambda]')
%     xlim([0 1])
%     ylabel('SD')
%     legend('Location','northeastoutside')
%     savefig(f1,[FigureFile1 '.fig'])
%     print(f1, '-dpsc','-fillpage', '-append', [FigureFile1 '.ps'])
%     
%     %% save figure2
%     figure(f2)
%     histogram(ax2,d_eq,bins,'DisplayStyle','stairs','Normalization','pdf','DisplayName',data_info.ID)
%     
%     
%     xlabel('distance [Lambda]')
%     xlim([0 1])
%     ylabel('SD')
%     legend('Location','northeastoutside')
%     savefig(f1,[FigureFile2 '.fig'])
%     print(f1, '-dpsc','-fillpage', '-append', [FigureFile2 '.ps'])
%     
%     %% save figure3
%     figure(f3)
%     histogram(ax3,d_eq,bins,'DisplayStyle','stairs','Normalization','pdf','DisplayName',data_info.ID)
%     
%     
%     xlabel('distance [Lambda]')
%     xlim([0 1])
%     ylabel('SD')
%     legend('Location','northeastoutside')
%     savefig(f1,[FigureFile3 '.fig'])
%     print(f1, '-dpsc','-fillpage', '-append', [FigureFile3 '.ps'])

    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     return 
end







