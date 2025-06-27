function PlotPwDensity(experiments,AnimalDataFolder,DataFolderMain,FigureFolder,PwDensityType,N_bootstrapSamples,BiasDataFile,colors)
    % %% parameter
    % experiment_Num = 9;
    % animal = 'dunnart';
    % PwDensityType ='WeightedPwDensityFixedFilter';
    % AnimalDataFolder='~/CIDBN/';
    % DataFolderMain = '~/Cloud/Cloud/PhD/MarsupialData/marsupial-data/DataTestGIF/';
    % FigureFolder = '~/Cloud/Cloud/PhD/Writing/phd_thesis/Dunnart_Paper/Figures/';

    %experiments.dunnart = [1 2 4 5 6 7 8 9];
    %PlotPwDensity(experiments,'~/Cloud/','~/Cloud/Cloud/PhD/MarsupialData/marsupial-data/DataHPC_GIF_adaptedFilter/')
    
    close all

    if nargin < 5
        FigureFolder = DataFolderMain;
    end
    if nargin < 6
        PwDensityType ='WeightedPwDensityFixedFilter';
    end
    if nargin < 7
        N_bootstrapSamples = 100;
    end

    if nargin < 8
        BiasDataFile = '';
        %BiasDataFolder = ['~/Cloud/Cloud/PhD/Writing/phd_thesis/OPM_Methods/MakeNoiseFromDataROI_ColumnSpacing/ResultDatav0/' animal '/'PwDensityBiasCorrection.mat''];
    end
    if ~isempty(BiasDataFile)
        load(BiasDataFile,'PwDensityBiasCorrection')
    end

    if nargin < 9
        colors = [ 255,0,255; 218,112,214]./255;
%         colors = [ 1 0 0; 0 0 1; 0 1 0; 1 0 1; 0 1 1; 0.5 0.5 0.5; 1 0.5 0; 0.5 0 0.5; 0 0.5 0.5];
    end

    MarkerSize = 6;
    LineWidth = 2;

    %% get animalList
    AnimalList = fieldnames(experiments);

    %% figure folder
    mkdir(FigureFolder)

    FigureFile = [FigureFolder PwDensityType];
    FigureFileSpacing = [FigureFolder 'ColumnSpacing'];


    %% collect colum spacing data for fit
    MarsupialList = {'wallaby'};
    PrimateList = {'mouselemur','galago','macaque','microcebus','shrew'};
    CarnivoreList = {'dunnart','cat','ferret','dunnart'};

    ColumnSpacing_Marsupial = [];
    BodyWeight_Marsupial = [];
    meanColumnSpacingMarssupials = [];
    meanBodyWeight_Marsupials = [];
    seColumnSpacingMarssupials = [];
    seBodyWeight_Marsupials = [];
    fit_weight_marsupial = [];

    ColumnSpacing_Carnivore = [];
    BodyWeight_Carnivore = [];
    meanColumnSpacingCarnivores = [];
    meanBodyWeight_Carnivores = [];
    seColumnSpacingCarnivores = [];
    seBodyWeight_Carnivores = [];
    fit_weight_carnivore = [];

    ColumnSpacing_Primate = [];
    BodyWeight_Primate = [];
    meanColumnSpacingPrimates = [];
    meanBodyWeight_Primates = [];
    seColumnSpacingPrimates = [];
    seBodyWeight_Primates = [];
    fit_weight_primates = [];
    
    %% make Figures
    f1 = figure();
    f1.Position = [1 1 800 300];
    ax1= axes('Parent',f1);


    f2 = figure();
    f2.Position = [1 1 800 300];
    ax2= axes('Parent',f2);


    for ii_animal = 1:length(AnimalList)

        animal = AnimalList{ii_animal};
        experiment_num_list = experiments.(animal);
    
        %% make figures animal
        f3 = figure();
        % f3.Position = [1 1 800 300];
        ax3= axes('Parent',f3);
    
        f4 = figure();
        ax4= axes('Parent',f4);

        f5 = figure(5);
        f5.Position = [1 1 300 300];
        ax5= axes('Parent',f5);
    
        %% pinwheel density sample
        MeanPwDensityBS = 0;
        PwDensityBS = zeros([max(experiment_num_list) N_bootstrapSamples])*NaN;
    %     PwDensityJS = zeros([max(experiment_num_list) size(data_obj.data,4)])*NaN;
        NumHypercolumnsBS = zeros([max(experiment_num_list) N_bootstrapSamples])*NaN;
        %SamplesJS = zeros([experiment_Num N_bootstrapSamples]);
        MeanPwDensities = [];
        NumHypercolumns = zeros([max(experiment_num_list) 1])*NaN;
        MeanPwDensities = zeros([max(experiment_num_list) 1])*NaN;
        weight_list = zeros([max(experiment_num_list) 1])*NaN;

        ColumnSpacing_Species = [];
        BodyWeight_Species = [];
    
    
        %% Plot Dunnart data
        for experiment_num = experiment_num_list
    
            %% data folder
            DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
    
            %% load animal info
            [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder);
    
            %% load spacing data
            SpacingFile = [DataFolder 'MapSpacingFiltered_' data_info.ID '.mat'];
            load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI')
            
            %% load CI spacing data
            CISpacingFile = [DataFolder 'CI_MapSpacingFiltered_' data_info.ID '.mat'];
            %load(CISpacingFile,'jackstat_average_spacing_mm','average_spacings_mm')
            load(CISpacingFile,'average_spacings_mm','local_spacings_mm','jackstat_average_spacing_mm','local_spacingsJS_mm','newROIsBS','newROIsJS')
            
            %% get spacing lacal mean, std and CIs
            newROI = (newROI==1);
            local_CS_mean = mean(local_spacing_mm(newROI));
            local_CS_std = std(local_spacing_mm(newROI));

            BS_local_CS_mean = zeros(sum(newROI(:)),N_bootstrapSamples-1);
            BS_local_CS_std = zeros(sum(newROI(:)),N_bootstrapSamples-1);
            for ii = 2:N_bootstrapSamples
                BS = local_spacings_mm{ii};
                BS_local_CS_mean(:,ii-1) = mean(BS(newROI));
                BS_local_CS_std(:,ii-1) = std(BS(newROI));
            end

            JS_local_CS_mean = zeros(sum(newROI(:)),size(local_spacingsJS_mm,2));
            JS_local_CS_std = zeros(sum(newROI(:)),size(local_spacingsJS_mm,2));
            for ii = 1:size(local_spacingsJS_mm,2)
                JS = local_spacingsJS_mm{ii};
                JS_local_CS_mean(:,ii) = mean(JS(newROI));
                JS_local_CS_std(:,ii) = std(JS(newROI));
            end

            CI_mean = bootstrap_ci(local_CS_mean,BS_local_CS_mean,jackstat_average_spacing_mm,1-0.68);
            CI_std = bootstrap_ci(local_CS_std,BS_local_CS_std,jackstat_average_spacing_mm,1-0.68);

            %% plot mean vs std with CIs for x and y axis in figure 5
            figure(f5)
            errorbar(ax5,local_CS_mean,local_CS_std, ...
                    local_CS_std-CI_std(1),CI_std(2)-local_CS_std, ...
                    local_CS_mean-CI_mean(1),CI_mean(2)-local_CS_mean,'ob')
            hold on;

    
            %% load pinwheel data
            PwInfoFile = [DataFolder 'PwInfo_' data_info.ID '.mat'];
            load(PwInfoFile,'PwInfo')
    
    
            %% load pinwheel CI data
            CIPwFile = [DataFolder 'CI_PwDensity_' data_info.ID '.mat'];
            load(CIPwFile,'alpha','PwInfosBS','PwInfosJS')
    
            %% fix naming for plotting
            data_info.ID = replace(data_info.ID,'_RightHemisphere',' RH');
    
            %% get pinwheel density
            MeanPwDensityCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,PwDensityType,false,alpha);
            if ~isempty(BiasDataFile)
                MeanPwDensity = [MeanPwDensityCI(1) PwInfo.(PwDensityType) MeanPwDensityCI(2)]-PwDensityBiasCorrection(experiment_Num);
                PwDensityBS(experiment_num,:) = getSamples(PwInfosBS,PwDensityType,data_obj,false)-PwDensityBiasCorrection(experiment_Num);
    %             PwDensityJS(experiment_num,:) = getSamples(PwInfosJS,PwDensityType,data_obj,false)-PwDensityBiasCorrection(experiment_Num);
            else
                MeanPwDensity = [MeanPwDensityCI(1) PwInfo.(PwDensityType) MeanPwDensityCI(2)];
                PwDensityBS(experiment_num,:) = getSamples(PwInfosBS,PwDensityType,data_obj,false);
    %             PwDensityJS(experiment_num,:) = getSamples(PwInfosJS,PwDensityType,data_obj,false);
            end
            disp([data_info.ID '  ' num2str(MeanPwDensity)])
            
    
            %% get mean spacing
            CI_average_spacing_mm = bootstrap_ci(average_spacings_mm(2:end),average_spacing_mm,jackstat_average_spacing_mm,alpha);
            MeanSpacing_mm = [CI_average_spacing_mm(1) average_spacing_mm CI_average_spacing_mm(2)];
            disp([num2str(MeanSpacing_mm) ' mm'])
    
            %% get animal infos
            weight_list(experiment_num) = data_obj.info.weight_in_grams;
            NumHypercolumns(experiment_num) = PwInfo.NumHypercolumns;
            MeanPwDensities(experiment_num) = MeanPwDensity(2); 
            NumHypercolumnsBS(experiment_num,:) = getSamples(PwInfosBS,'NumHypercolumns',data_obj,false);
    
    
            %% plot mean pinwheel density
            figure(f1)
            plot(ax1,[data_obj.info.weight_in_grams],[PwInfo.(PwDensityType)],'o','DisplayName',data_info.ID,'Color',colors(ii_animal,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
            hold on
            
    
            %% plot mean spacing
            figure(f2)
            plot(ax2,[data_obj.info.weight_in_grams],[average_spacing_mm],'o','DisplayName',data_info.ID,'Color',colors(ii_animal,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
            hold on
    
            %% plot pw density CI
            figure(f3)
            errorbar(ax3,experiment_num,[MeanPwDensity(2)], ...
                    [MeanPwDensity(2)-MeanPwDensity(1)], ...
                    [MeanPwDensity(3)-MeanPwDensity(2)],'o','DisplayName',data_info.ID)
            hold on
    
    
            %% plot column spacing CI
            figure(f4)
            errorbar(ax4,experiment_num,[MeanSpacing_mm(2)], ...
                    [MeanSpacing_mm(2)-MeanSpacing_mm(1)], ...
                    [MeanSpacing_mm(3)-MeanSpacing_mm(2)],'o','DisplayName',data_info.ID)
            hold on

            %% add data to list for fit
            if any(strcmp(CarnivoreList,animal))
                ColumnSpacing_Carnivore = [ColumnSpacing_Carnivore; average_spacing_mm];
                BodyWeight_Carnivore = [BodyWeight_Carnivore; data_obj.info.weight_in_grams];
                fit_weight_carnivore = [fit_weight_carnivore; 1/length(experiment_num_list)];
            elseif any(strcmp(PrimateList,animal))
                ColumnSpacing_Primate = [ColumnSpacing_Primate; average_spacing_mm];
                BodyWeight_Primate = [BodyWeight_Primate; data_obj.info.weight_in_grams];
                fit_weight_primates = [fit_weight_primates; 1/length(experiment_num_list)];
            elseif any(strcmp(MarsupialList,animal))
                ColumnSpacing_Marsupial = [ColumnSpacing_Marsupial; average_spacing_mm];
                BodyWeight_Marsupial = [BodyWeight_Marsupial; data_obj.info.weight_in_grams];
                fit_weight_marsupial = [fit_weight_marsupial; 1/length(experiment_num_list)];
            end

            ColumnSpacing_Species = [ColumnSpacing_Species; average_spacing_mm];
            BodyWeight_Species = [BodyWeight_Species; data_obj.info.weight_in_grams];
    
            
        end

        if any(strcmp(CarnivoreList,animal))
            meanColumnSpacingCarnivores = [meanColumnSpacingCarnivores mean(ColumnSpacing_Carnivore)];
            meanBodyWeight_Carnivores = [meanBodyWeight_Carnivores mean(BodyWeight_Carnivore)];
            seColumnSpacingCarnivores = [seColumnSpacingCarnivores std(ColumnSpacing_Carnivore)/sqrt(length(ColumnSpacing_Carnivore))];
            seBodyWeight_Carnivores = [seBodyWeight_Carnivores std(BodyWeight_Carnivore)/sqrt(length(BodyWeight_Carnivore))];
        elseif any(strcmp(PrimateList,animal))  
            meanColumnSpacingPrimates = [meanColumnSpacingPrimates mean(ColumnSpacing_Primate)];
            meanBodyWeight_Primates = [meanBodyWeight_Primates mean(BodyWeight_Primate)];
            seColumnSpacingPrimates = [seColumnSpacingPrimates std(ColumnSpacing_Primate)/sqrt(length(ColumnSpacing_Primate))];
            seBodyWeight_Primates = [seBodyWeight_Primates std(BodyWeight_Primate)/sqrt(length(BodyWeight_Primate))];
        elseif any(strcmp(MarsupialList,animal))
            meanColumnSpacingMarssupials = [meanColumnSpacingMarssupials mean(ColumnSpacing_Marsupial)];
            meanBodyWeight_Marsupials = [meanBodyWeight_Marsupials mean(BodyWeight_Marsupial)];
            seColumnSpacingMarssupials = [seColumnSpacingMarssupials std(ColumnSpacing_Marsupial)/sqrt(length(ColumnSpacing_Marsupial))];
            seBodyWeight_Marsupials = [seBodyWeight_Marsupials std(BodyWeight_Marsupial)/sqrt(length(BodyWeight_Marsupial))];
        end



        %% add mean with CI to fig1
        figure(1)
        plot_bca_line(MeanPwDensities, weight_list,'Color',colors(ii_animal,:),'Alpha', 1-0.68,'nbootstrap',1000)
        hold on

        %% Finish Plots and save figures
        LineWidth=3;

    
        %% calc mean corrected pinwheel density BS all animals
        figure(f3)
        % MeanPwDensityAnimals = MeanPwDensityAnimals/experiment_Num;
        %PwDensityCI = bootstrap_ci(SamplesBS,MeanPwDensityAnimals,SamplesJS,alpha);
        % PwDensityCI(1)=prctile(SamplesBS,2.5);
        % PwDensityCI(2)=prctile(SamplesBS,100-2.5);
        % disp(std(SamplesBS)/sqrt(length(SamplesBS))*1.96)
        PwDensityBS=PwDensityBS(~isnan(PwDensityBS)); 
        SortedPwDensityBS = sort(PwDensityBS);
        MeanPwDensityBS = mean(PwDensityBS,'all');
        PwDensityCI(1)= SortedPwDensityBS(floor(length(SortedPwDensityBS)*alpha/2));
        PwDensityCI(2)= SortedPwDensityBS(ceil(length(SortedPwDensityBS)*(1-alpha/2)));
        errorbar(ax3,0,[MeanPwDensityBS], ...
                [PwDensityCI(1)-MeanPwDensityBS], ...
                [PwDensityCI(2)-MeanPwDensityBS],'or','DisplayName','average')
        disp([PwDensityCI(1) MeanPwDensityBS PwDensityCI(2)])
        
    %     %% calc mean corrected pinwheel density BS animals
    %     MeanPwDensityAnimals = mean(MeanPwDensities,'all');
    %     samplesBS = BootstrapSamplesMean1D(MeanPwDensities,N_bootstrapSamples,1);
    %     samplesJK = JackknifeSamplesMean1D(MeanPwDensities);
    %     PwDensityCI = bootstrap_ci(samplesBS,mean(MeanPwDensities),samplesJK,alpha);
    %     errorbar(ax3,0,[MeanPwDensityAnimals], ...
    %             [PwDensityCI(1)-MeanPwDensityAnimals], ...
    %             [PwDensityCI(2)-MeanPwDensityAnimals],'xr')
    %     disp([PwDensityCI(1) MeanPwDensityAnimals PwDensityCI(2)])
    
    
        %% calc mean with NumHypercolumns weight
        NumHypercolumns = NumHypercolumns(~isnan(MeanPwDensities));
        MeanPwDensities = MeanPwDensities(~isnan(MeanPwDensities));
        weight = NumHypercolumns./mean(NumHypercolumns,'all');
        
        MeanPwDensities=MeanPwDensities.*weight;
        
        
        SortedPwDensity = sort(MeanPwDensities);
        MeanPwDensity = mean(MeanPwDensities,'all');
    %     PwDensityCI(1)= SortedPwDensity(flo(length(SortedPwDensity)*alpha/2));
    %     PwDensityCI(2)= SortedPwDensity(ceil(length(SortedPwDensity)*(1-alpha/2)));
    %     errorbar(ax3,0,[MeanPwDensity], ...
    %             [PwDensityCI(1)-MeanPwDensity], ...
    %             [PwDensityCI(2)-MeanPwDensity],'+r')
    %     disp([PwDensityCI(1) MeanPwDensityBS PwDensityCI(2)])
        disp([MeanPwDensity std(SortedPwDensity)])


    
        %% save figure3
        figure(f3)
        hold on
        plot([0 max(experiment_num_list)+1], [pi pi],"DisplayName","\pi")
        title('Pinwheel Density')
        xlim([0 max(experiment_num_list)+1])
        set(gca,'DefaultLineLineWidth',LineWidth)
        ylim([0 8])
        ylabel('ρ [1/Λ²]')
        xlabel('animal')
        legend('Location','northeastoutside')
    
    
        savefig(f3,[FigureFile 'CI.fig'])
        print(f3,'-depsc2', [FigureFile animal 'CI.eps']);
    
    
        %% save figure4
        figure(f4)
        hold on
        title('Column Spacing')
        xlim([0 max(experiment_num_list)+1])
        set(gca,'DefaultLineLineWidth',LineWidth)
        ylim([0 1.5])
        ylabel('Λ [mm]')
        xlabel('animal')
        legend('Location','northeastoutside')
    
        savefig(f4,[FigureFileSpacing 'CI.fig'])
        print(f4,'-depsc2', [FigureFileSpacing animal 'CI.eps']);


        %% save figure5
        figure(f5)
        hold on
        title('Local Column Spacing Mean and standard deviation')
        set(gca,'DefaultLineLineWidth',LineWidth)
        xlabel('mean Λ [mm]')
        ylabel('std Λ [mm]')

        savefig(f5,[FigureFileSpacing 'LocalSpacingMeanVar.fig'])
        print(f5,'-depsc2', [FigureFileSpacing animal 'LocalSpacingMeanVar.eps']);
    
    
        %% plot corr Number of Hypercolumns vs MeanPwDensity
        figure
        plot(NumHypercolumns,MeanPwDensities,'*')
        xlabel('Number Hypercolumns')
        ylabel('mean pw density')
        print('-depsc2', [FigureFile animal 'PwDensVsHypecolumnN.eps']);
    end



    %% plot data wallaby paper
    %  - **Wallaby B** - 2.37, Λ of 1.74 mm, **Weight: 3900 g**
    % - **Wallaby C** - 1.95, Λ of 1.14 mm, **Weight: 3700 g**
    % - **Wallaby D** - 2.95, Λ of 1.27 mm, **Weight: 4000 g**
    % - **Wallaby E** - 2.86, Λ of 1.61 mm, **Weight: 4200 g**
    % - **Wallaby F** - 2.97, Λ of 1.55 mm, **Weight: 4200 g**
    % - **Wallaby H** - 2.65, Λ of 1.28 mm, **Weight: 5700 g**

    wallaby_weight_g = [3900 3700 4000 4200 4200 5700];
    wallaby_pw_density = [2.37 1.95 2.95 2.86 2.97 2.65];
    wallaby_colum_spacing_mm = [1.74 1.14 1.27 1.61 1.55 1.28];


    % pw density
    figure(f1)
    plot(ax1,wallaby_weight_g,wallaby_pw_density,'o','DisplayName','Wallaby','Color',colors(length(AnimalList)+1,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
    hold on
    plot_bca_line(wallaby_pw_density, wallaby_weight_g,'Color',colors(length(AnimalList)+1,:),'Alpha', 1-0.68,'nbootstrap',1000)
    hold on

    % column spacing
    figure(f2)
    plot(ax2,wallaby_weight_g,wallaby_colum_spacing_mm,'o','DisplayName','Wallaby','Color',colors(length(AnimalList)+1,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
    hold on

    % add data to list for fit
    ColumnSpacing_Marsupial = [ColumnSpacing_Marsupial; wallaby_colum_spacing_mm'];
    BodyWeight_Marsupial = [BodyWeight_Marsupial; wallaby_weight_g'];
    fit_weight_marsupial = [fit_weight_marsupial; 1/length(wallaby_weight_g)*ones(length(wallaby_weight_g),1)];

    meanColumnSpacingMarssupials = [meanColumnSpacingMarssupials mean(wallaby_colum_spacing_mm)];
    meanBodyWeight_Marsupials = [meanBodyWeight_Marsupials mean(wallaby_weight_g)];
    seColumnSpacingMarssupials = [seColumnSpacingMarssupials std(wallaby_colum_spacing_mm)/sqrt(length(wallaby_colum_spacing_mm))];
    seBodyWeight_Marsupials = [seBodyWeight_Marsupials std(wallaby_weight_g)/sqrt(length(wallaby_weight_g))];

    %% Load and plot data Mouse lemur paper
    file_path = '/home/michael/Cloud/Cloud/PhD/Writing/phd_thesis/Dunnart_Paper/DataMouseLemurPaper/';
    addpath(file_path)
    DataMouseLemurPaper = read_data_tables(file_path);

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.ferret,[255,210,131],MarkerSize,"o",LineWidth);%colors(length(AnimalList)+2,:)
    ColumnSpacing_Carnivore = [ColumnSpacing_Carnivore; DataMouseLemurPaper.ferret.WL_mm];
    BodyWeight_Carnivore = [BodyWeight_Carnivore; DataMouseLemurPaper.ferret.Weight_g];
    fit_weight_carnivore = [fit_weight_carnivore; 1/length(DataMouseLemurPaper.ferret.Weight_g)*ones(length(DataMouseLemurPaper.ferret.Weight_g),1)];

    Weight_g = DataMouseLemurPaper.ferret.Weight_g;
    WL_mm = DataMouseLemurPaper.ferret.WL_mm;
    idx_NoNan = ~isnan(Weight_g)&~isnan(WL_mm);
    Weight_g = Weight_g(idx_NoNan);
    WL_mm = WL_mm(idx_NoNan);
    meanColumnSpacingCarnivores = [meanColumnSpacingCarnivores mean(WL_mm)];
    meanBodyWeight_Carnivores = [meanBodyWeight_Carnivores mean(Weight_g)];
    seColumnSpacingCarnivores = [seColumnSpacingCarnivores std(WL_mm)/sqrt(length(WL_mm))];
    seBodyWeight_Carnivores = [seBodyWeight_Carnivores std(Weight_g)/sqrt(length(Weight_g))];
   

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.cat,[230,220,229],MarkerSize,"o",LineWidth);%colors(length(AnimalList)+3,:)
    ColumnSpacing_Carnivore = [ColumnSpacing_Carnivore; DataMouseLemurPaper.cat.WL_mm];
    BodyWeight_Carnivore = [BodyWeight_Carnivore; DataMouseLemurPaper.cat.Weight_g];
    fit_weight_carnivore = [fit_weight_carnivore; 1/length(DataMouseLemurPaper.cat.Weight_g)*ones(length(DataMouseLemurPaper.cat.Weight_g),1)];

    Weight_g = DataMouseLemurPaper.cat.Weight_g;
    WL_mm = DataMouseLemurPaper.cat.WL_mm;
    idx_NoNan = ~isnan(Weight_g)&~isnan(WL_mm);
    Weight_g = Weight_g(idx_NoNan);
    WL_mm = WL_mm(idx_NoNan);
    meanColumnSpacingCarnivores = [meanColumnSpacingCarnivores mean(WL_mm)];
    meanBodyWeight_Carnivores = [meanBodyWeight_Carnivores mean(Weight_g)];
    seColumnSpacingCarnivores = [seColumnSpacingCarnivores std(WL_mm)/sqrt(length(WL_mm))];
    seBodyWeight_Carnivores = [seBodyWeight_Carnivores std(Weight_g)/sqrt(length(Weight_g))];

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.macaque,[0,0,255],MarkerSize,"d",LineWidth);
    ColumnSpacing_Primate = [ColumnSpacing_Primate; DataMouseLemurPaper.macaque.WL_mm];
    BodyWeight_Primate = [BodyWeight_Primate; DataMouseLemurPaper.macaque.Weight_g];
    fit_weight_primates = [fit_weight_primates; 1/length(DataMouseLemurPaper.macaque.Weight_g)*ones(length(DataMouseLemurPaper.macaque.Weight_g),1)];

    Weight_g = DataMouseLemurPaper.macaque.Weight_g;
    WL_mm = DataMouseLemurPaper.macaque.WL_mm;
    idx_NoNan = ~isnan(Weight_g)&~isnan(WL_mm);
    Weight_g = Weight_g(idx_NoNan);
    WL_mm = WL_mm(idx_NoNan);
    meanColumnSpacingPrimates = [meanColumnSpacingPrimates mean(WL_mm)];
    meanBodyWeight_Primates = [meanBodyWeight_Primates mean(Weight_g)];
    seColumnSpacingPrimates = [seColumnSpacingPrimates std(WL_mm)/sqrt(length(WL_mm))];
    seBodyWeight_Primates = [seBodyWeight_Primates std(Weight_g)/sqrt(length(Weight_g))];

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.microcebus,[255,0,0],MarkerSize,"d",LineWidth);
    ColumnSpacing_Primate = [ColumnSpacing_Primate; DataMouseLemurPaper.microcebus.WL_mm];
    BodyWeight_Primate = [BodyWeight_Primate; DataMouseLemurPaper.microcebus.Weight_g];
    fit_weight_primates = [fit_weight_primates; 1/length(DataMouseLemurPaper.microcebus.Weight_g)*ones(length(DataMouseLemurPaper.microcebus.Weight_g),1)];

    Weight_g = DataMouseLemurPaper.microcebus.Weight_g;
    WL_mm = DataMouseLemurPaper.microcebus.WL_mm;
    idx_NoNan = ~isnan(Weight_g)&~isnan(WL_mm);
    Weight_g = Weight_g(idx_NoNan);
    WL_mm = WL_mm(idx_NoNan);   
    meanColumnSpacingPrimates = [meanColumnSpacingPrimates mean(WL_mm)];
    meanBodyWeight_Primates = [meanBodyWeight_Primates mean(Weight_g)];
    seColumnSpacingPrimates = [seColumnSpacingPrimates std(WL_mm)/sqrt(length(WL_mm))];
    seBodyWeight_Primates = [seBodyWeight_Primates std(Weight_g)/sqrt(length(Weight_g))];

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.galago,[50,139,160],MarkerSize,"d",LineWidth);
    ColumnSpacing_Primate = [ColumnSpacing_Primate; DataMouseLemurPaper.galago.WL_mm];
    BodyWeight_Primate = [BodyWeight_Primate; DataMouseLemurPaper.galago.Weight_g];
    fit_weight_primates = [fit_weight_primates; 1/length(DataMouseLemurPaper.galago.Weight_g)*ones(length(DataMouseLemurPaper.galago.Weight_g),1)];

    Weight_g = DataMouseLemurPaper.galago.Weight_g;
    WL_mm = DataMouseLemurPaper.galago.WL_mm;
    idx_NoNan = ~isnan(Weight_g)&~isnan(WL_mm);
    Weight_g = Weight_g(idx_NoNan);
    WL_mm = WL_mm(idx_NoNan);
    meanColumnSpacingPrimates = [meanColumnSpacingPrimates mean(WL_mm)];
    meanBodyWeight_Primates = [meanBodyWeight_Primates mean(Weight_g)]; 
    seColumnSpacingPrimates = [seColumnSpacingPrimates std(WL_mm)/sqrt(length(WL_mm))];
    seBodyWeight_Primates = [seBodyWeight_Primates std(Weight_g)/sqrt(length(Weight_g))];

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.shrew,[69,178,135],MarkerSize,"d",LineWidth);
    ColumnSpacing_Primate = [ColumnSpacing_Primate; DataMouseLemurPaper.shrew.WL_mm];
    BodyWeight_Primate = [BodyWeight_Primate; DataMouseLemurPaper.shrew.Weight_g];
    fit_weight_primates = [fit_weight_primates; 1/length(DataMouseLemurPaper.shrew.Weight_g)*ones(length(DataMouseLemurPaper.shrew.Weight_g),1)];

    Weight_g = DataMouseLemurPaper.shrew.Weight_g;
    WL_mm = DataMouseLemurPaper.shrew.WL_mm;
    idx_NoNan = ~isnan(Weight_g)&~isnan(WL_mm);
    Weight_g = Weight_g(idx_NoNan);
    WL_mm = WL_mm(idx_NoNan);
    meanColumnSpacingPrimates = [meanColumnSpacingPrimates mean(WL_mm)];
    meanBodyWeight_Primates = [meanBodyWeight_Primates mean(Weight_g)];
    seColumnSpacingPrimates = [seColumnSpacingPrimates std(WL_mm)/sqrt(length(WL_mm))];
    seBodyWeight_Primates = [seBodyWeight_Primates std(Weight_g)/sqrt(length(Weight_g))];

    %% save figure1
    figure(f1)
    WeightRange = [1 11*10^3];%[(min(weight_list)*0.9) (max(weight_list)*1.1)];
    plot(ax1,WeightRange,(WeightRange-WeightRange+1)*pi,'--black','DisplayName','\pi')

    title(PwDensityType)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'DefaultLineLineWidth',LineWidth)
    xlim([10 11*10^3])
    xticks(10.^(0:4))
    ylim([0 6])
    ylabel('Pinwheel Density')
    xlabel('Body Weight [g]')
    legend('Location','northeastoutside')


    savefig(f1,[FigureFile '.fig'])
    %print(f1,'-depsc2', [FigureFile,'.eps']);
    print(f1, '-dpdf', [FigureFile, '.pdf']);

    %% add power-law fits and save figure2
    figure(f2)

    % remove nans from fit data
    idx_noNan = ~isnan(ColumnSpacing_Carnivore)&~isnan(BodyWeight_Carnivore);
    ColumnSpacing_Carnivore_woNan = ColumnSpacing_Carnivore(idx_noNan);
    BodyWeight_Carnivore_woNan = BodyWeight_Carnivore(idx_noNan);
    fit_weight_carnivore_woNan = fit_weight_carnivore(idx_noNan);

    idx_noNan = ~isnan(ColumnSpacing_Primate)&~isnan(BodyWeight_Primate);
    ColumnSpacing_Primate_woNan = ColumnSpacing_Primate(idx_noNan);
    BodyWeight_Primate_woNan = BodyWeight_Primate(idx_noNan);
    fit_weight_primates_woNan = fit_weight_primates(idx_noNan);

    idx_noNan = ~isnan(ColumnSpacing_Marsupial)&~isnan(BodyWeight_Marsupial);
    ColumnSpacing_Marsupial_woNan = ColumnSpacing_Marsupial(idx_noNan);
    BodyWeight_Marsupial_woNan = BodyWeight_Marsupial(idx_noNan);
    fit_weight_marsupial_woNan = fit_weight_marsupial(idx_noNan);


    % fit parameter
    nboot = 1000;
    % confidence equivalent to SE
    conf_level = 68;
    x_range = 10.^(0:0.01:4);
    x_range = reshape(x_range,[length(x_range) 1]);
% 
    % fit carnivore with BcA fit_powerlaw_bca and plot
    [mean_params_carnivore, y_fit_car, y_lo_car, y_hi_car] = fit_logarithmic_bca_curveci(BodyWeight_Carnivore_woNan, ColumnSpacing_Carnivore_woNan, fit_weight_carnivore_woNan, x_range, conf_level, nboot);
    plot_logarithmic_curveci(x_range, y_fit_car, y_lo_car, y_hi_car, 'LineStyle','-')
    hold on
    annotate_logarithmic_fit_point(mean_params_carnivore, 25, 'mouse weasel')


    % fit primate with BcA fit_powerlaw_bca and plot
    [mean_params, y_fit, y_lo, y_hi] = fit_logarithmic_bca_curveci(BodyWeight_Primate_woNan, ColumnSpacing_Primate_woNan, fit_weight_primates_woNan, x_range, conf_level, nboot);
    plot_logarithmic_curveci(x_range, y_fit, y_lo, y_hi, 'LineStyle','--')

    %plot_powerlaw_curveci_noc(x_eval, y_fit, y_lo, y_hi, varargin)

%     %  fit carnivore
%     [xFitCarnivore,yFitCarnivore,fitresultCarnivore] = fitPowerLaw(BodyWeight_Carnivore_woNan,ColumnSpacing_Carnivore_woNan,100,[10 10000]);
%     plot(ax2,xFitCarnivore,yFitCarnivore,'--black','DisplayName','Carnivore Fit')
% 
%     % fit primate
%     [xFitPrimate,yFitPrimate,fitresultPrimate] = fitPowerLaw(BodyWeight_Primate_woNan,ColumnSpacing_Primate_woNan,100,[10 10000]);
%     plot(ax2,xFitPrimate,yFitPrimate,':black','DisplayName','Primate Fit')

    % % fit marsupial
    % [xFitMarsupial,yFitMarsupial,fitresultMarsupial] = fitPowerLaw(BodyWeight_Marsupial_woNan,ColumnSpacing_Marsupial_woNan);
    % plot(ax2,xFitMarsupial,yFitMarsupial,'-.black','DisplayName','Marsupial Fit')

    %title('Column Spacing [mm]')
    set(gca,'DefaultLineLineWidth',LineWidth)
    set(gca,'xscale','log')
    %set(gca,'yscale','log')
    xlim([10 11*10^3])
    xticks(10.^(0:4))
    ylim([0 1.8])
    ylabel('Column Spacing [mm]')
    xlabel('Body Weight [g]')
    legend('Location','northeastoutside')

    %% inset
    inset_ax = axes('Position', [0.15 0.65 0.25 0.25]); 
    box on
    plot_logarithmic_curveci(x_range, y_fit_car, y_lo_car, y_hi_car, 'LineStyle','-')
    hold on
    annotate_logarithmic_fit_point(mean_params_carnivore, 4.2, 'long-tailed planigale')
    hold on
    annotate_logarithmic_fit_point(mean_params_carnivore, 2, 'bumblebee bat')
    xlim([1 10])
    set(inset_ax, 'FontSize', 8) % Optional: make inset font smaller
    set(inset_ax,'xscale','log')
    %xticks(1:10)



    savefig(f2,[FigureFileSpacing '.fig'])
    %print(f2,'-depsc2', [FigureFileSpacing,'.eps']);
    print(f2, '-dpdf', [FigureFileSpacing, '.pdf']);

    %% add marsupials to carnivores
    meanColumnSpacingCarnivores = [meanColumnSpacingCarnivores meanColumnSpacingMarssupials];
    meanBodyWeight_Carnivores = [meanBodyWeight_Carnivores meanBodyWeight_Marsupials];
    seColumnSpacingCarnivores = [seColumnSpacingCarnivores seColumnSpacingMarssupials];
    seBodyWeight_Carnivores = [seBodyWeight_Carnivores seBodyWeight_Marsupials];



    %% get max and min weights from mean lists
    minWeight = min([meanBodyWeight_Primates meanBodyWeight_Carnivores])*0.9;
    maxWeight = max([meanBodyWeight_Primates meanBodyWeight_Carnivores])*1.1;
    %% save data for fit mean
    figure;

    % plot primate data with 2d error bars
    errorbar(meanBodyWeight_Primates,meanColumnSpacingPrimates,seColumnSpacingPrimates,seColumnSpacingPrimates,seBodyWeight_Primates,seBodyWeight_Primates,'o','DisplayName','Primate data')
    hold on

    % plot carnivore data with 2d error bars
    errorbar(meanBodyWeight_Carnivores,meanColumnSpacingCarnivores,seColumnSpacingCarnivores,seColumnSpacingCarnivores,seBodyWeight_Carnivores,seBodyWeight_Carnivores,'o','DisplayName','Carnivore & Marsupial data')
    hold on


     % fit and plot primate
    %[xFitPrimate,yFitPrimate,fitresultPrimate] = weightedPowerlawFit(meanBodyWeight_Primates,meanColumnSpacingPrimates,seBodyWeight_Primates,seColumnSpacingPrimates,10.^(0:0.1:4));
    [xFitPrimate,yFitPrimate,fitresultPrimate] = fitPowerLaw(meanBodyWeight_Primates,meanColumnSpacingPrimates,100,[minWeight maxWeight]);
    plot(xFitPrimate,yFitPrimate,':black','DisplayName','Primate Fit')
    hold on

    % fit and plot carnivore
    %[xFitCarnivore,yFitCarnivore,fitresultCarnivore] = weightedPowerlawFit(meanBodyWeight_Carnivores,meanColumnSpacingCarnivores,seBodyWeight_Carnivores,seColumnSpacingCarnivores,10.^(0:0.1:4));
    [xFitCarnivore,yFitCarnivore,fitresultCarnivore] = fitPowerLaw(meanBodyWeight_Carnivores,meanColumnSpacingCarnivores,100,[minWeight maxWeight]);
    plot(xFitCarnivore,yFitCarnivore,'--black','DisplayName','Carnivore Fit')
    hold on

    % set x-axes to log scale
    set(gca,'xscale','log')

    % set labels
    xlim([1 11*10^3])
    xticks(10.^(0:4))
    ylim([0 1.8])
    ylabel('Column Spacing [mm]')
    xlabel('Body Weight [g]')
    legend()

    savefig([FigureFileSpacing 'Mean.fig'])
    print('-depsc2', [FigureFileSpacing 'AnimalMeanColumnSpacing.eps']);


    
    
    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

function [x_fit,y_fit,p_fit] = weightedPowerlawFit(x, y, sigma_x, sigma_y, x_fit)
    % POWERLAW_FIT Performs a power-law fit with uncertainties in both x and y.
    %
    % Inputs:
    %   x - Independent variable values
    %   y - Dependent variable values
    %   sigma_x - Uncertainties in x
    %   sigma_y - Uncertainties in y
    %   x_fit - Values of x for which to compute the fitted curve
    %
    % Outputs:
    %   A_fit - Fitted coefficient A
    %   B_fit - Fitted exponent B
    %   y_fit - Computed fitted values for x_fit
    %
    % Define power-law function
    powerlaw = @(p, x) p(1) * x.^p(2); % y = A*x^B
    
    % Initial guess for parameters A and B
    p0 = [1, 1];
    
    % Define the objective function (sum of squared residuals)
    obj_fun = @(p) ((powerlaw(p, x) - y) ./ sigma_y).^2 + ((x - (y / p(1)).^(1 / p(2))) ./ sigma_x).^2;
    
    % Perform nonlinear least squares fitting
    options = optimset('Display', 'off', 'TolFun', 1e-8);
    p_fit = lsqnonlin(@(p) sqrt(obj_fun(p)), p0, [], [], options);
    
    % Extract fitted parameters
    A_fit = p_fit(1);
    B_fit = p_fit(2);
    
    % Compute fitted curve for given x_fit values
    y_fit = A_fit * x_fit.^B_fit;
    
    % Display results
    fprintf('Power-law fit: y = %.4f * x^{%.4f}\n', A_fit, B_fit);
end

function [xFit,yFit,fitresult] = fitPowerLaw(x,y,nPoints,xlim)

    if nargin < 3
        nPoints = 100;
    end
    if nargin<4
        xlim=[min(x) max(x)];
    end

    % Define the power-law fit model: y = a*x^b
    powerLawModel = fittype('a*x^b', 'independent', 'x', 'coefficients', {'a', 'b'});

    % Perform the fit
    fitresult = fit(x(:), y(:), powerLawModel);

    % % Display fit parameters
    % disp(fitresult)

    % Generate smooth curve for plotting
    xFit = linspace(xlim(1),xlim(2), nPoints); % More points for a smooth curve
    yFit = fitresult.a * xFit.^fitresult.b;
end

function [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,data,color,MarkerSize,MarkerSymbol,LineWidth)

    color = color./255;

    data.PwD

    figure(f1)
    plot(ax1,data.Weight_g,data.PwD,'DisplayName',data.legend_text,'Color',color, 'Marker', MarkerSymbol, 'MarkerSize', MarkerSize, 'LineStyle', 'none', 'LineWidth', LineWidth)
    hold on

    plot_bca_line(data.PwD, data.Weight_g,'Color',color,'Alpha', 1-0.68,'nbootstrap',1000)
    hold on
   

    figure(f2)
    plot(ax2,data.Weight_g,data.WL_mm,'DisplayName',data.legend_text,'Color',color, 'Marker', MarkerSymbol, 'MarkerSize', MarkerSize, 'LineStyle', 'none', 'LineWidth', LineWidth)
    hold on

end

function [mean_params, ci_bca] = fit_powerlaw_bca(x, y, w, conf_level, nboot)
% Fits y = a*x^b with weights and computes BCa bootstrap confidence intervals
% Inputs:
%   x, y        - Data vectors
%   w           - Weights
%   conf_level  - Confidence level (e.g., 95)
%   nboot       - Number of bootstrap samples
% Outputs:
%   mean_params - [a, b] from fit to original data
%   ci_bca      - [a_lo a_hi; b_lo b_hi] (BCa intervals)

    x = x(:); y = y(:); w = w(:);
    n = numel(x);
    alpha = (100 - conf_level) / 2 / 100;

    % Define model (no constant term)
    model = @(b, x) b(1) * x.^b(2);  % [a, b]
    beta0 = [range(y)/range(x.^1), 1];  % reasonable initial guess

    % Initial weighted fit
    sqrtw = sqrt(w);
    modelw = @(b,x) model(b,x) .* sqrtw;
    [theta_hat, ~, ~, ~, ~] = nlinfit(x, y .* sqrtw, modelw, beta0);

    %% Bootstrap samples
    boot_params = NaN(nboot, 2);
    for i = 1:nboot
        idx = randsample(n, n, true);
        xb = x(idx); yb = y(idx); wb = w(idx);
        sqrtwb = sqrt(wb);
        modelwb = @(b,x) model(b,x) .* sqrtwb;
        try
            boot_params(i,:) = nlinfit(xb, yb .* sqrtwb, modelwb, beta0);
        catch
            % leave NaN
        end
    end
    boot_params = boot_params(all(~isnan(boot_params),2), :);
    mean_params = theta_hat;

    %% Jackknife samples
    jack_params = NaN(n, 2);
    for i = 1:n
        idx = setdiff(1:n, i);
        xj = x(idx); yj = y(idx); wj = w(idx);
        sqrtwj = sqrt(wj);
        modelwj = @(b,x) model(b,x) .* sqrtwj;
        try
            jack_params(i,:) = nlinfit(xj, yj .* sqrtwj, modelwj, beta0);
        catch
            % leave NaN
        end
    end
    jack_params = jack_params(all(~isnan(jack_params),2), :);

    %% Compute BCa confidence intervals
    ci_bca = NaN(2, 2);
    for j = 1:2  % loop over parameters: a, b
        boot_theta = boot_params(:,j);
        jack_theta = jack_params(:,j);
        theta0 = theta_hat(j);

        z0 = norminv(sum(boot_theta < theta0) / length(boot_theta));
        jack_mean = mean(jack_theta);
        diffs = jack_mean - jack_theta;
        skew = sum(diffs.^3) / (6 * (sum(diffs.^2).^(3/2)));
        a = skew;

        z_low = norminv(alpha);
        z_high = norminv(1 - alpha);
        pct_low = normcdf(z0 + (z0 + z_low) / (1 - a * (z0 + z_low))) * 100;
        pct_high = normcdf(z0 + (z0 + z_high) / (1 - a * (z0 + z_high))) * 100;

        ci = prctile(boot_theta, [pct_low, pct_high]);
        ci_bca(j,:) = ci;
    end
end



function plot_powerlaw_bca_fit(mean_params, ci_bca, x_range, varargin)
% Plots y = a*x^b with BCa confidence band over specified x_range
% Inputs:
%   mean_params - [a, b] from fit_powerlaw_bca_noc
%   ci_bca      - [a_lo a_hi; b_lo b_hi], BCa confidence intervals
%   x_range     - Vector of x values to evaluate and plot
%
% Optional (name-value pairs):
%   'Alpha'     - Transparency of CI band (default: 0.2)
%   'Color'     - Color for fit and CI (default: black)
%   'LineStyle' - Line style for fit line (default: '-')
%   'Data'      - [x_data, y_data] to overlay original data

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'Alpha', 0.2);
    addParameter(p, 'Color', [0 0 0]);
    addParameter(p, 'LineStyle', '-');
    addParameter(p, 'Data', []);
    parse(p, varargin{:});
    alpha = p.Results.Alpha;
    col = p.Results.Color;
    linestyle = p.Results.LineStyle;
    data = p.Results.Data;

    % Extract fit parameters
    a = mean_params(1);
    b = mean_params(2);

    % Evaluate fit
    x_fit = x_range(:);
    y_fit = a * x_fit.^b;

    % Evaluate lower and upper CI curves
    a_lo = ci_bca(1,1); a_hi = ci_bca(1,2);
    b_lo = ci_bca(2,1); b_hi = ci_bca(2,2);

    y_lo = a_lo * x_fit.^b_lo;
    y_hi = a_hi * x_fit.^b_hi;

    % Plot shaded CI area
    fill([x_fit; flipud(x_fit)], [y_lo; flipud(y_hi)], col, ...
         'FaceAlpha', alpha, 'EdgeColor', 'none'); hold on;

    % Plot mean fit line
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2, 'LineStyle', linestyle);

    % Optional data overlay
    if ~isempty(data)
        scatter(data(:,1), data(:,2), 25, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    end
end

function [mean_params, y_fit, y_lo, y_hi] = fit_powerlaw_bca_curveci_noc(x, y, w, x_eval, conf_level, nboot)
% Fits y = a*x^b (no constant) and computes BCa bootstrap confidence intervals on the curve
% Inputs:
%   x, y, w     - Data and weights
%   x_eval      - Points to evaluate the fit and CIs
%   conf_level  - Confidence level (e.g., 95)
%   nboot       - Number of bootstrap samples
% Outputs:
%   mean_params - [a, b] from original fit
%   y_fit       - Fit values at x_eval
%   y_lo, y_hi  - Lower and upper BCa bounds at x_eval

    x = x(:); y = y(:); w = w(:); x_eval = x_eval(:);
    n = numel(x);
    alpha = (100 - conf_level) / 2;

    % Model: y = a * x^b
    model = @(b, x) b(1) * x.^b(2);
    beta0 = [range(y)/range(x.^1), 1];

    % Fit original model
    sqrtw = sqrt(w);
    modelw = @(b,x) model(b,x) .* sqrtw;
    [theta_hat, ~, ~, ~, ~] = nlinfit(x, y .* sqrtw, modelw, beta0);
    mean_params = theta_hat;
    y_fit = model(mean_params, x_eval);

    %% Bootstrap predictions
    y_boot = NaN(nboot, numel(x_eval));
    for i = 1:nboot
        idx = randsample(n, n, true);
        xb = x(idx); yb = y(idx); wb = w(idx);
        sqrtwb = sqrt(wb);
        modelwb = @(b,x) model(b,x) .* sqrtwb;
        try
            b_fit = nlinfit(xb, yb .* sqrtwb, modelwb, beta0);
            y_boot(i,:) = model(b_fit, x_eval);
        catch
            % leave NaN
        end
    end
    y_boot = y_boot(all(~isnan(y_boot), 2), :);
    nboot_valid = size(y_boot, 1);

    %% Jackknife predictions
    y_jack = NaN(n, numel(x_eval));
    for i = 1:n
        idx = setdiff(1:n, i);
        xj = x(idx); yj = y(idx); wj = w(idx);
        sqrtwj = sqrt(wj);
        modelwj = @(b,x) model(b,x) .* sqrtwj;
        try
            b_jack = nlinfit(xj, yj .* sqrtwj, modelwj, beta0);
            y_jack(i,:) = model(b_jack, x_eval);
        catch
            % leave NaN
        end
    end
    valid_jack = all(~isnan(y_jack), 2);
    y_jack = y_jack(valid_jack, :);

    %% BCa interval per x_eval
    y_lo = NaN(size(x_eval));
    y_hi = NaN(size(x_eval));
    for k = 1:numel(x_eval)
        yb = y_boot(:,k);
        yj = y_jack(:,k);
        y0 = y_fit(k);

        % Bias correction
        z0 = norminv(sum(yb < y0) / nboot_valid);

        % Acceleration
        jack_mean = mean(yj);
        diffs = jack_mean - yj;
        skew = sum(diffs.^3) / (6 * (sum(diffs.^2).^(3/2)));
        a = skew;

        % Adjusted percentiles
        z_low = norminv(alpha / 100);
        z_high = norminv(1 - alpha / 100);
        pct_low = normcdf(z0 + (z0 + z_low) / (1 - a * (z0 + z_low))) * 100;
        pct_high = normcdf(z0 + (z0 + z_high) / (1 - a * (z0 + z_high))) * 100;

        ci = prctile(yb, [pct_low, pct_high]);
        y_lo(k) = ci(1);
        y_hi(k) = ci(2);
    end
end

function plot_powerlaw_curveci_noc(x_eval, y_fit, y_lo, y_hi, varargin)
% Plots y = a*x^b with BCa confidence bands
    p = inputParser;
    addParameter(p, 'Alpha', 0.2);
    addParameter(p, 'Color', [0 0 0]);
    addParameter(p, 'LineStyle', '-');
    addParameter(p, 'Data', []);
    parse(p, varargin{:});
    alpha = p.Results.Alpha;
    col = p.Results.Color;
    linestyle = p.Results.LineStyle;
    data = p.Results.Data;

    % Plot CI band
    fill([x_eval; flipud(x_eval)], [y_lo(:); flipud(y_hi(:))], col, ...
         'FaceAlpha', alpha, 'EdgeColor', 'none'); hold on;

    % Plot mean fit
    plot(x_eval, y_fit, 'Color', col, 'LineWidth', 2, 'LineStyle', linestyle);

    % Overlay data if provided
    if ~isempty(data)
        scatter(data(:,1), data(:,2), 25, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    end
end


function [mean_params, y_fit, y_lo, y_hi] = fit_logarithmic_bca_curveci(x, y, w, x_eval, conf_level, nboot)
% Fits y = a*log(x) + b with weights and computes BCa CIs of the curve
% Inputs:
%   x, y, w     - Data vectors and weights
%   x_eval      - Evaluation points for CI curves
%   conf_level  - Confidence level (e.g., 95)
%   nboot       - Number of bootstrap samples
% Outputs:
%   mean_params - [a, b] from original fit
%   y_fit       - Mean fit at x_eval
%   y_lo, y_hi  - BCa confidence bounds at x_eval

    x = x(:); y = y(:); w = w(:); x_eval = x_eval(:);
    n = numel(x);
    alpha = (100 - conf_level) / 2;

    % Define model
    model = @(b, x) b(1) * log(x) + b(2);  % y = a*log(x) + b
    beta0 = [1, mean(y)];

    % Fit to original data
    sqrtw = sqrt(w);
    modelw = @(b,x) model(b,x) .* sqrtw;
    [theta_hat, ~, ~, ~, ~] = nlinfit(x, y .* sqrtw, modelw, beta0);
    mean_params = theta_hat;
    y_fit = model(mean_params, x_eval);

    %% Bootstrap predictions
    y_boot = NaN(nboot, numel(x_eval));
    for i = 1:nboot
        idx = randsample(n, n, true);
        xb = x(idx); yb = y(idx); wb = w(idx);
        sqrtwb = sqrt(wb);
        modelwb = @(b,x) model(b,x) .* sqrtwb;
        try
            b_fit = nlinfit(xb, yb .* sqrtwb, modelwb, beta0);
            y_boot(i,:) = model(b_fit, x_eval);
        catch
            % leave NaN row
        end
    end
    y_boot = y_boot(all(~isnan(y_boot), 2), :);
    nboot_valid = size(y_boot, 1);

    %% Jackknife predictions
    y_jack = NaN(n, numel(x_eval));
    for i = 1:n
        idx = setdiff(1:n, i);
        xj = x(idx); yj = y(idx); wj = w(idx);
        sqrtwj = sqrt(wj);
        modelwj = @(b,x) model(b,x) .* sqrtwj;
        try
            b_jack = nlinfit(xj, yj .* sqrtwj, modelwj, beta0);
            y_jack(i,:) = model(b_jack, x_eval);
        catch
            % leave NaN
        end
    end
    valid_jack = all(~isnan(y_jack), 2);
    y_jack = y_jack(valid_jack, :);

    %% BCa interval per x_eval
    y_lo = NaN(size(x_eval));
    y_hi = NaN(size(x_eval));
    for k = 1:numel(x_eval)
        yb = y_boot(:,k);
        yj = y_jack(:,k);
        y0 = y_fit(k);

        % Bias correction
        z0 = norminv(sum(yb < y0) / nboot_valid);

        % Acceleration
        jack_mean = mean(yj);
        diffs = jack_mean - yj;
        skew = sum(diffs.^3) / (6 * (sum(diffs.^2).^(3/2)));
        a = skew;

        % Adjusted percentiles
        z_low = norminv(alpha / 100);
        z_high = norminv(1 - alpha / 100);
        pct_low = normcdf(z0 + (z0 + z_low) / (1 - a * (z0 + z_low))) * 100;
        pct_high = normcdf(z0 + (z0 + z_high) / (1 - a * (z0 + z_high))) * 100;

        ci = prctile(yb, [pct_low, pct_high]);
        y_lo(k) = ci(1);
        y_hi(k) = ci(2);
    end
end

function plot_logarithmic_curveci(x_eval, y_fit, y_lo, y_hi, varargin)
% Plots y = a*log(x) + b with BCa confidence bands
% Inputs:
%   x_eval - x values used for evaluation
%   y_fit  - predicted values (mean fit)
%   y_lo   - lower BCa bound
%   y_hi   - upper BCa bound

    p = inputParser;
    addParameter(p, 'Alpha', 0.2);
    addParameter(p, 'Color', [0 0 0]);
    addParameter(p, 'LineStyle', '-');
    addParameter(p, 'Data', []);
    parse(p, varargin{:});
    alpha = p.Results.Alpha;
    col = p.Results.Color;
    linestyle = p.Results.LineStyle;
    data = p.Results.Data;

    % Plot CI band
    fill([x_eval; flipud(x_eval)], [y_lo(:); flipud(y_hi(:))], col, ...
         'FaceAlpha', alpha, 'EdgeColor', 'none'); hold on;

    % Plot mean fit
    plot(x_eval, y_fit, 'Color', col, 'LineWidth', 2, 'LineStyle', linestyle);

    % Overlay data if present
    if ~isempty(data)
        scatter(data(:,1), data(:,2), 25, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    end
end


function plot_bca_line(data, x_values, varargin)
% Plot the mean and BCa confidence interval as a horizontal line with shaded area.
%
% Inputs:
%   data      - Vector of data points
%   x_values  - Vector of x positions (length ≥ 2). Will span [min(x_values), max(x_values)].
%
% Optional name-value pairs:
%   'Alpha'        - Significance level (default: 0.05)
%   'NBootstrap'   - Number of bootstrap samples (default: 10000)
%   'StatFun'      - Statistic function (default: @mean)
%   'Color'        - Line and fill color (default: [0 0 0])
%   'LineWidth'    - Line width (default: 2)
%   'FaceAlpha'    - Shaded area transparency (default: 0.2)

    % Default values
    alpha = 0.05;
    n_bootstrap = 10000;
    stat_fun = @mean;
    line_color = [0 0 0];
    line_width = 2;
    face_alpha = 0.2;

    % Parse optional arguments
    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'alpha'
                alpha = varargin{i+1};
            case 'nbootstrap'
                n_bootstrap = varargin{i+1};
            case 'statfun'
                stat_fun = varargin{i+1};
            case 'color'
                line_color = varargin{i+1};
            case 'linewidth'
                line_width = varargin{i+1};
            case 'facealpha'
                face_alpha = varargin{i+1};
        end
    end

    % Ensure x_values is a row vector and get min/max
    x_values = x_values(:)';
    x_min = min(x_values);
    x_max = max(x_values);

    % Compute mean and BCa CI
    [stat_val, ci_lower, ci_upper] = bca_confidence_interval(data, alpha, n_bootstrap, stat_fun);

    % Plot shaded CI area
    hold on;
    fill([x_min, x_max, x_max, x_min], ...
         [ci_lower, ci_lower, ci_upper, ci_upper], ...
         line_color, 'FaceAlpha', face_alpha, 'EdgeColor', 'none');

    % Plot mean line
    plot([x_min, x_max], [stat_val, stat_val], '-', ...
         'Color', line_color, 'LineWidth', line_width);
end

function [stat_original, ci_lower, ci_upper] = bca_confidence_interval(data, alpha, n_bootstrap, stat_fun)
% Compute the BCa confidence interval for a statistic (default: mean)
% using bootstrap and jackknife resampling.
%
% Inputs:
%   data         - Vector of input data
%   alpha        - Significance level (e.g., 0.05 for 95% CI)
%   n_bootstrap  - Number of bootstrap samples
%   stat_fun     - Function handle for the statistic (default: @mean)
%
% Outputs:
%   stat_original - Original statistic
%   ci_lower      - Lower bound of BCa confidence interval
%   ci_upper      - Upper bound of BCa confidence interval

    if nargin < 4
        stat_fun = @mean;
    end

    data = data(:);
    n = length(data);
    stat_original = stat_fun(data);

    % Bootstrap resampling
    boot_stats = zeros(n_bootstrap, 1);
    for i = 1:n_bootstrap
        sample = data(randi(n, n, 1));
        boot_stats(i) = stat_fun(sample);
    end
    boot_stats = sort(boot_stats);

    % Bias-correction
    z0 = norminv(mean(boot_stats < stat_original));

    % Jackknife resampling
    jack_stats = zeros(n, 1);
    for i = 1:n
        jack_sample = data([1:i-1, i+1:end]);
        jack_stats(i) = stat_fun(jack_sample);
    end
    jack_mean = mean(jack_stats);
    jack_dev = jack_mean - jack_stats;

    % Acceleration
    num = sum(jack_dev.^3);
    denom = 6 * (sum(jack_dev.^2))^(3/2);
    a = num / denom;

    % Adjusted confidence levels
    z_alpha_low = norminv(alpha / 2);
    z_alpha_high = norminv(1 - alpha / 2);

    pct_low = normcdf(z0 + (z0 + z_alpha_low) / (1 - a * (z0 + z_alpha_low)));
    pct_high = normcdf(z0 + (z0 + z_alpha_high) / (1 - a * (z0 + z_alpha_high)));

    % Convert percentiles to indexes
    idx_low = max(1, floor(pct_low * n_bootstrap));
    idx_high = min(n_bootstrap, ceil(pct_high * n_bootstrap));

    ci_lower = boot_stats(idx_low);
    ci_upper = boot_stats(idx_high);
end

function annotate_logarithmic_fit_point(params, x_val, label, varargin)
% Annotates a point on the logarithmic fit y = a*log(x) + b
% Inputs:
%   params - [a, b] from the logarithmic fit
%   x_val  - Scalar x value to evaluate the fit at
%   label  - String label to display next to the point
%
% Optional name-value pairs:
%   'Color'  - Marker and text color (default: 'k')
%   'Offset' - [dx, dy] offset for the text label in data units (default: [0.05, 0.05])

    p = inputParser;
    addParameter(p, 'Color', 'k');
    addParameter(p, 'Offset', [0.05, 0.05]);
    parse(p, varargin{:});

    col = p.Results.Color;
    offset = p.Results.Offset;

    % Compute fit value
    a = params(1);
    b = params(2);
    y_val = a * log(x_val) + b;

    % Plot marker at the point
    plot(x_val, y_val, 'o', 'Color', col, ...
         'MarkerFaceColor', col, 'MarkerSize', 6); hold on;

    % Plot text label slightly offset
    text(x_val + offset(1), y_val + offset(2), label, ...
         'Color', col, 'FontSize', 10, 'Interpreter', 'none');
end



