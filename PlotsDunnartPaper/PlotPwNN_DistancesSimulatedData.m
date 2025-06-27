function PlotPwNN_DistancesSimulatedData(animal_list,experiment_Num,ResultDataMain,FigureFolder)
    % experiment_Num.dunnart = [1 2 5 6 7];
    % PlotPwNN_DistancesSimulatedData(["dunnart"],experiment_Num,'~/Cloud/Cloud/PhD/Writing/phd_thesis/OPM_Methods/MakeNoiseFromDataROI_ColumnSpacing/ResultDatav1kTestROI11Mar25/')
    



    %% parameter
    experiment_Num = checkFormatNum(experiment_Num);
    if nargin <4
        FigureFolder = ResultDataMain;
    end


    
    %% figure folder
    mkdir(FigureFolder)

    
    for ii = 1:length(animal_list)
        animal = animal_list(ii);
        animal_=convertStringsToChars(animal);
        

        ResultData = [ResultDataMain animal_ '/'];

        %% Distance Lists
        d = [];
        d_eq = [];
        d_op= [];
        all_areas = [];
        all_n = [];

        d_ROI = [];
        d_eq_ROI = [];
        d_op_ROI= [];
        all_areas_ROI = [];
        all_n_ROI = [];


        if ~isa(experiment_Num,'struct') && (length(experiment_Num) == length(animal_list))
                experiment_Num = 1: experiment_Num(ii);
        else
            experiment_Num_list = experiment_Num.(animal_);
        end
        for experiment_num = experiment_Num_list

            
            load([ResultData num2str(experiment_num) 'PwDataROI_ColumnSpacing.mat'],...
                'ColumnSpacin_pix_woROI_woNoise','d_woROI_woNoise','d_eq_woROI_woNoise','d_op_woROI_woNoise','circ_areas_woROI_woNoise','nPws_woROI_woNoise',...
                'ColumnSpacin_pix_ROI_woNoise','d_ROI_woNoise','d_eq_ROI_woNoise','d_op_ROI_woNoise','circ_areas_ROI_woNoise','nPws_ROI_woNoise');

            %% flatten  array
            ColumnSpacin_pix_woROI_woNoise = ColumnSpacin_pix_woROI_woNoise(experiment_num,:);
            ColumnSpacin_pix_ROI_woNoise = ColumnSpacin_pix_ROI_woNoise(experiment_num,:);

            %% flatten cell arrays to one big array
            d_woROI_woNoise = cell2mat(d_woROI_woNoise(experiment_num,:));
            d_eq_woROI_woNoise = cell2mat(d_eq_woROI_woNoise(experiment_num,:));
            d_op_woROI_woNoise = cell2mat(d_op_woROI_woNoise(experiment_num,:));
            circ_areas_woROI_woNoise = cell2mat(circ_areas_woROI_woNoise(experiment_num,:));
            nPws_woROI_woNoise = cell2mat(nPws_woROI_woNoise(experiment_num,:));

            d_ROI_woNoise = cell2mat(d_ROI_woNoise(experiment_num,:));
            d_eq_ROI_woNoise = cell2mat(d_eq_ROI_woNoise(experiment_num,:));
            d_op_ROI_woNoise = cell2mat(d_op_ROI_woNoise(experiment_num,:));
            circ_areas_ROI_woNoise = cell2mat(circ_areas_ROI_woNoise(experiment_num,:));
            nPws_ROI_woNoise = cell2mat(nPws_ROI_woNoise(experiment_num,:));

            
            %% caclac NN distances
            d = [d d_woROI_woNoise/mean(ColumnSpacin_pix_woROI_woNoise)];
            d_eq = [d_eq d_eq_woROI_woNoise/mean(ColumnSpacin_pix_woROI_woNoise)];
            d_op = [d_op d_op_woROI_woNoise/mean(ColumnSpacin_pix_woROI_woNoise)];
            all_areas = [all_areas circ_areas_woROI_woNoise];
            all_n = [all_n nPws_woROI_woNoise];

            d_ROI = [d_ROI d_ROI_woNoise/mean(ColumnSpacin_pix_ROI_woNoise)];
            d_eq_ROI = [d_eq_ROI d_eq_ROI_woNoise/mean(ColumnSpacin_pix_ROI_woNoise)];
            d_op_ROI = [d_op_ROI d_op_ROI_woNoise/mean(ColumnSpacin_pix_ROI_woNoise)];
            all_areas_ROI = [all_areas_ROI circ_areas_ROI_woNoise];
            all_n_ROI = [all_n_ROI nPws_ROI_woNoise];

        end

        distances.(animal_).d = d;
        distances.(animal_).d_eq = d_eq;
        distances.(animal_).d_op= d_op;
        distances.(animal_).all_areas = all_areas;
        distances.(animal_).all_n = all_n;

        distances.(animal_).d_ROI = d_ROI;
        distances.(animal_).d_eq_ROI = d_eq_ROI;
        distances.(animal_).d_op_ROI= d_op_ROI;
        distances.(animal_).all_areas_ROI = all_areas_ROI;
        distances.(animal_).all_n_ROI = all_n_ROI;
        

    end


    %% load comparison data
%     microcebus = load('microcebus_Huber.mat');
%     macaque = load('macaque_Angelucci.mat');
    
    %% plot NN distances
    edges = 0:0.05:1;
    f = figure();
    f.Position = [100 100 2000 500];
    t = tiledlayout(1,4);
    
    %% d
    nexttile; %figure(1);
    plot_nn_theory('d')
    for ii = 1:length(animal_list)
        animal_ = char(animal_list(ii));
        hold on
        plot_nn_data_dist(distances.(animal_).d,edges, animal_)
        hold on
        plot_nn_data_dist(distances.(animal_).d_ROI,edges, [animal_ ' ROI'])
    end
% 
%     if ~isempty(MainFolderNoiseSimulation)
%         for ii = 1:length(animal_list)
%             animal_ = animal_list(ii);
%             hold on
%             plot_nn_data_dist(distances.(animal_).d_noise,edges, [char(animal) ' Noise'])
%         end
%     end
%     hold on
%     plot_nn_data_dist(microcebus.d,edges, 'Microcebus')
%     hold on
%     plot_nn_data_dist(macaque.d,edges, 'Macaque')
    
    title('d')
    xlim([0 1])
    ylim([0 1])
    xlabel('NN distance [\Lambda]')
    ylabel('Freq. (normalized)')
    set(gca,'xtick',[0 0.5 1])
    set(gca,'ytick',[0 0.25 0.5 0.75 1])
    set(gca,'yticklabel',{'0','','0.5','','1'})
    box off
    set(gca,'Fontsize',25)
    set(gca,'linewidth',1.5)
    legend()
    
    %print('-depsc2', [FigureFolder 'comparison_nndist.eps']);
    
    
    %% d**
    nexttile;%figure(2);
    plot_nn_theory('d++')
    for ii = 1:length(animal_list)
        animal_ = char(animal_list(ii));
        hold on
        plot_nn_data_dist(distances.(animal_).d_eq,edges, animal_)
        hold on
        plot_nn_data_dist(distances.(animal_).d_eq_ROI,edges, [animal_ ' ROI'])
    end
    
    % if ~isempty(MainFolderNoiseSimulation)
    %     for ii = 1:length(animal_list)
    %         animal_ = animal_list(ii);
    %         hold on
    %         plot_nn_data_dist(distances.(animal_).d_eq_noise,edges, [char(animal) ' Noise'])
    %     end
    % end

    %hold on
    %plot_nn_data_dist(d_eq,edges, 'Dunnart')
%     hold on
%     plot_nn_data_dist(microcebus.d_eq,edges, 'Microcebus')
%     hold on
%     plot_nn_data_dist(macaque.d_eq,edges, 'Macaque')
    
    title('d++')
    xlim([0 1])
    ylim([0 1])
    xlabel('NN distance [\Lambda]')
    ylabel('Freq. (normalized)')
    set(gca,'xtick',[0 0.5 1])
    set(gca,'ytick',[0 0.25 0.5 0.75 1])
    set(gca,'yticklabel',{'0','','0.5','','1'})
    box off
    set(gca,'Fontsize',25)
    set(gca,'linewidth',1.5)
    
    %print('-depsc2', [FigureFolder 'comparison_nndist_eq.eps']);
    
    
    %% d**
    nexttile;%figure(3);
    plot_nn_theory('d+-')
    for ii = 1:length(animal_list)
        animal_ = char(animal_list(ii));
        hold on
        plot_nn_data_dist(distances.(animal_).d_op,edges, animal_)
        hold on
        plot_nn_data_dist(distances.(animal_).d_op_ROI,edges, [animal_ ' ROI'])
    end

    % if ~isempty(MainFolderNoiseSimulation)
    %     for ii = 1:length(animal_list)
    %         animal_ = animal_list(ii);
    %         hold on
    %         plot_nn_data_dist(distances.(animal_).d_op_noise,edges, [char(animal) ' Noise'])
    %     end
    % end

%     hold on
%     plot_nn_data_dist(microcebus.d_op,edges, 'Microcebus')
%     hold on
%     plot_nn_data_dist(macaque.d_op,edges, 'Macaque')
    
    title('d+-')
    xlim([0 1])
    ylim([0 1])
    xlabel('NN distance [\Lambda]')
    ylabel('Freq. (normalized)')
    set(gca,'xtick',[0 0.5 1])
    set(gca,'ytick',[0 0.25 0.5 0.75 1])
    set(gca,'yticklabel',{'0','','0.5','','1'})
    box off
    set(gca,'Fontsize',25)
    set(gca,'linewidth',1.5)
    
    %print('-depsc2', [FigureFolder 'comparison_nndist_op.eps']);
    
    %% plot count sd
    nexttile;
    plot_count_sd_theory()
    for ii = 1:length(animal_list)
        animal_ = char(animal_list(ii));
        hold on
        plot_count_sd_data(distances.(animal_).all_areas, distances.(animal_).all_n, animal_)
        hold on
        plot_count_sd_data(distances.(animal_).all_areas_ROI, distances.(animal_).all_n_ROI, [animal_ ' ROI'])
    end

%     hold on
%     plot_count_sd_data(microcebus.all_areas, microcebus.all_n, 'Microcebus')
%     hold on
%     plot_count_sd_data(macaque.all_areas, macaque.all_n, 'Macaque')
    
    title('SD')
    xlim([0.3 20])
    ylim([0.1 5])
    ylabel('SD of pinwheel counts')
    xlabel('Region size [\Lambda]')
    set(gca,'xtick',[0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20])
    set(gca,'xticklabel',{'','','','','','','','1.0','','','','','','','','','10',''})
    set(gca,'ytick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5])
    set(gca,'yticklabel',{'0.1','','','','','','','','','1.0','','','',''})
    set(gca,'Fontsize',25)
    
    
    %% save figure
%     if BS
%         print(f,'-depsc2', [FigureFolder 'comparison_nndistBS.eps']);
%     else
%         
% 
%     end
    print(f,'-depsc2', [FigureFolder 'TheoryDataComparison_nndist.eps']);
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







