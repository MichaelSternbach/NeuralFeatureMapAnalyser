function PlotPwNN_DistancesMultiAnimal(animal_list,experiment_Num,AnimalDataFolder,DataFolderMain,BS,FigureFolder,MainFolderNoiseSimulation)
    % experiment_Num.dunnart = [2 3 4 5 6 7 9];
    % PlotPwNN_DistancesMultiAnimal(["dunnart"],experiment_Num,'~/CIDBN/','~/Cloud/Cloud/PhD/MarsupialData/marsupial-data/DataHPC_GIF/')
    % '~/Cloud/Cloud/PhD/Writing/phd_thesis/OPM_Methods/Figures/','~/Cloud/Cloud/PhD/Writing/phd_thesis/OPM_Methods/MakeNoiseFromDataROI_ColumnSpacing/ResultDatav1k/')
    



    %% parameter
    experiment_Num = checkFormatNum(experiment_Num);
    if nargin <5
        BS = false;
    end
    if nargin <6
        FigureFolder = DataFolderMain;
    end
    if nargin <7
        MainFolderNoiseSimulation = '';
    end

%     PwDensityType= 'MeanPwDensity';
    
    %% figure folder
    mkdir(FigureFolder)

    
    for ii = 1:length(animal_list)
        animal = animal_list(ii);
        animal_=convertStringsToChars(animal);
        %% Distance Lists
        d = [];
        d_eq = [];
        d_op= [];
        all_areas = [];
        all_n = [];

        if ~isa(experiment_Num,'struct') && (length(experiment_Num) == length(animal_list))
                experiment_Num = 1: experiment_Num(ii);
        else
            experiment_Num_list = experiment_Num.(animal_);
        end
        for experiment_num = experiment_Num_list

            %% data folder
            DataFolder = [DataFolderMain lower(animal_) '/' lower(animal_) num2str(experiment_num) '/'];

            %% animal data 
            [data_info,~,data_obj,~,~] = getAnimalData(animal_,experiment_num,AnimalDataFolder);

            %% get column spacing
            %[average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,DataFolder,false);
            SpacingFile = [DataFolder 'MapSpacingFiltered_' data_obj.info.ID '.mat'];
            load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI','WavletCoefficient')
            

            %% load pinwheel data

            PwInfoFile = [DataFolder 'PwInfo_' data_obj.info.ID '.mat'];
            load(PwInfoFile,'PwInfo')

%             getCI = false;
%             do_plotting=0;
%             llp_cutoffs = linspace(0.01, 1,100);
%             beta=0.5;
%             Bootstrapsamples = 100;
%             data_obj.prepare_samples_array(Bootstrapsamples)
%             PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta);
            
            %% calc scaling
            scale = average_spacing_mm*data_info.pix_per_mm;%
            
            %% caclac NN distances
            d = [d PwInfo.d./scale];
            d_eq = [d_eq PwInfo.d_eq./scale];
            d_op = [d_op PwInfo.d_op./scale];


            %% pw area SDV
            all_areas = [all_areas PwInfo.circ_areas./data_info.pix_per_mm^2]; %^2 ./scale^2
            all_n = [all_n PwInfo.n];

            if BS
                %% load CI spacing data
                CISpacingFile = [DataFolder 'CI_MapSpacingFiltered_' data_obj.info.ID '.mat'];
                load(CISpacingFile,'average_spacings_mm')
                
        
                %% load pinwheel CI data
                CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
                load(CIPwFile,'alpha','PwInfosBS','PwInfosJS')
                
                %% loop over BS
                for ii_BS = 1: length(PwInfosBS)
%                     PwInfosBS{ii_BS}

                    %% calc scaling
                    scale = average_spacings_mm(ii_BS)*data_info.pix_per_mm;%
                    
                    try
                        %% caclac NN distances
                        d = [d PwInfosBS{ii_BS}.d./scale];
                        d_eq = [d_eq PwInfosBS{ii_BS}.d_eq./scale];
                        d_op = [d_op PwInfosBS{ii_BS}.d_op./scale];
                    catch
                        disp([animal num2str(experiment_num) ' BS' num2str(ii_BS) ' failed for NN calculation!'])
                    end

        
                    try
                        %% pw area SDV
                        all_areas = [all_areas PwInfosBS{ii_BS}.circ_areas./data_info.pix_per_mm^2]; %^2 ./scale^2
                        all_n = [all_n PwInfosBS{ii_BS}.n];
                    catch
                        disp([animal num2str(experiment_num) ' BS' num2str(ii_BS) ' failed for SVD calculation!'])
                    end

                end

            end



        end

        distances.(animal_).d = d;
        distances.(animal_).d_eq = d_eq;
        distances.(animal_).d_op= d_op;
        distances.(animal_).all_areas = all_areas;
        distances.(animal_).all_n = all_n;

        %% get NoiseSimulation data
        d = [];
        d_eq = [];
        d_op= [];
        all_areas = [];
        all_n = [];

        if ~isempty(MainFolderNoiseSimulation)
            NoiseSimulationFolder = [MainFolderNoiseSimulation lower(animal_) '/'];
            for experiment_num = experiment_Num_list
                NoiseSimulationFile = [NoiseSimulationFolder num2str(ii) 'PwDataROI_ColumnSpacing.mat'];
                load(NoiseSimulationFile,'d_ROI_Noise','d_eq_ROI_Noise','d_op_ROI_Noise')
                scale = average_spacing_mm*data_info.pix_per_mm;%
                
                d_ROI_Noise = flattenCellArray(d_ROI_Noise);
                d = [d d_ROI_Noise./scale];
                d_eq_ROI_Noise = flattenCellArray(d_eq_ROI_Noise);
                d_eq = [d_eq d_eq_ROI_Noise./scale];
                d_op_ROI_Noise = flattenCellArray(d_op_ROI_Noise);
                d_op = [d_op d_op_ROI_Noise./scale];
            end
            distances.(animal_).d_noise = d;
            distances.(animal_).d_eq_noise = d_eq;
            distances.(animal_).d_op_noise = d_op;

        end
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
        animal_ = animal_list(ii);
        hold on
        plot_nn_data_dist(distances.(animal_).d,edges, animal_)
    end

    if ~isempty(MainFolderNoiseSimulation)
        for ii = 1:length(animal_list)
            animal_ = animal_list(ii);
            hold on
            plot_nn_data_dist(distances.(animal_).d_noise,edges, [char(animal) ' Noise'])
        end
    end
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
        animal_ = animal_list(ii);
        hold on
        plot_nn_data_dist(distances.(animal_).d_eq,edges, animal_)
    end
    
    if ~isempty(MainFolderNoiseSimulation)
        for ii = 1:length(animal_list)
            animal_ = animal_list(ii);
            hold on
            plot_nn_data_dist(distances.(animal_).d_eq_noise,edges, [char(animal) ' Noise'])
        end
    end

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
        animal_ = animal_list(ii);
        hold on
        plot_nn_data_dist(distances.(animal_).d_op,edges, animal_)
    end

    if ~isempty(MainFolderNoiseSimulation)
        for ii = 1:length(animal_list)
            animal_ = animal_list(ii);
            hold on
            plot_nn_data_dist(distances.(animal_).d_op_noise,edges, [char(animal) ' Noise'])
        end
    end

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
        animal_ = animal_list(ii);
        hold on
        plot_count_sd_data(distances.(animal_).all_areas, distances.(animal_).all_n, animal_)
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
    if BS
        print(f,'-depsc2', [FigureFolder 'comparison_nndistBS.eps']);
    else
        print(f,'-depsc2', [FigureFolder 'comparison_nndist.eps']);

    end
    
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







