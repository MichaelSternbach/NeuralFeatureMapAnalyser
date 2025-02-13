function PlotPwDensity(experiment_num_list,animal,AnimalDataFolder,DataFolderMain,FigureFolder,PwDensityType,N_bootstrapSamples,BiasDataFile,colorlist)
    % %% parameter
    % experiment_Num = 9;
    % animal = 'dunnart';
    % PwDensityType ='WeightedPwDensityFixedFilter';
    % AnimalDataFolder='~/CIDBN/';
    % DataFolderMain = '~/Cloud/Cloud/PhD/MarsupialData/marsupial-data/DataTestGIF/';
    % FigureFolder = '~/Cloud/Cloud/PhD/Writing/phd_thesis/Dunnart_Paper/Figures/';

    if nargin < 5
        FigureFolder = [DataFolderMain animal '/'];
    end
    if nargin < 6
        PwDensityType ='WeightedPwDensityFixedFilter';;
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
        colorlist = ["r","g","b","c","m","y"];;
    end

    %% figure folder
    mkdir(FigureFolder)

    FigureFile = [FigureFolder PwDensityType];
    FigureFileSpacing = [FigureFolder 'ColumnSpacing'];

    %     rm_cmd = ['rm -f ' FigureFile '.ps'];
    %     disp(rm_cmd)
    %     system(rm_cmd)


    %% make Figures
    f1 = figure();
    f1.Position = [1 1 800 300];
    ax1= axes('Parent',f1);



    f2 = figure();
    f2.Position = [1 1 800 300];
    ax2= axes('Parent',f2);
    weight_list = [];



    f3 = figure();
    % f3.Position = [1 1 800 300];
    ax3= axes('Parent',f3);

    f4 = figure();
    ax4= axes('Parent',f4);

    %% pinwheel density sample
    MeanPwDensityAnimals = 0;
    SamplesBS = zeros([max(experiment_num_list) N_bootstrapSamples])*NaN;
    %SamplesJS = zeros([experiment_Num N_bootstrapSamples]);
    MeanPwDensities = [];
    NumHypercolumns = zeros([max(experiment_num_list) 1])*NaN;
    MeanPwDensities = zeros([max(experiment_num_list) 1])*NaN;
    weight_list = zeros([max(experiment_num_list) 1])*NaN;


    %% Plot Dunnart data
    ii=1;
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
        load(CISpacingFile,'jackstat_average_spacing_mm','average_spacings_mm')


        %% load pinwheel data
        PwInfoFile = [DataFolder 'PwInfo_' data_info.ID '.mat'];
        load(PwInfoFile,'PwInfo')


        %% load pinwheel CI data
        CIPwFile = [DataFolder 'CI_PwDensity_' data_info.ID '.mat'];
        load(CIPwFile,'alpha','PwInfosBS','PwInfosJS')

        %% fix naming for plotting
        data_info.ID = replace(data_info.ID,'_',' ');

        %% get pinwheel density
        MeanPwDensityCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,PwDensityType,false,alpha);
        if ~isempty(BiasDataFile)
            MeanPwDensity = [MeanPwDensityCI(1) PwInfo.(PwDensityType) MeanPwDensityCI(2)]-PwDensityBiasCorrection(experiment_Num);
        else
            MeanPwDensity = [MeanPwDensityCI(1) PwInfo.(PwDensityType) MeanPwDensityCI(2)];
        end
        disp([data_info.ID '  ' num2str(MeanPwDensity)])

        % %% get sample
        % BS = getSamples(PwInfosBS,PwDensityType,data_obj,false)-PwDensityBiasCorrection(experiment_Num);
        % MeanPwDensityAnimals = MeanPwDensityAnimals+BS(1);
        % %MeanPwDensities = [MeanPwDensities PwInfo.(PwDensityType)-PwDensityBiasCorrection(experiment_Num)];
        % SamplesBS(experiment_num,:)= BS;
        

        %% get mean spacing
        CI_average_spacing_mm = bootstrap_ci(average_spacings_mm(2:end),average_spacing_mm,jackstat_average_spacing_mm,alpha);
        MeanSpacing_mm = [CI_average_spacing_mm(1) average_spacing_mm CI_average_spacing_mm(2)];
        disp([num2str(MeanSpacing_mm) ' mm'])

        %% get animal infos
        weight_list(experiment_num) = data_obj.info.weight_in_grams;
        NumHypercolumns(experiment_num) = PwInfo.NumHypercolumns;
        MeanPwDensities(experiment_num) = MeanPwDensity(2); 

        %% plot mean pinwheel density
        figure(f1)
        plot(ax1,[data_obj.info.weight_in_grams],[PwInfo.(PwDensityType)],'o','DisplayName',data_info.ID,'Color',colorlist(ii))
        hold on

        %% plot mean spacing
        figure(f2)
        plot(ax2,[data_obj.info.weight_in_grams],[average_spacing_mm],'o','DisplayName',data_info.ID,'Color',colorlist(ii))
        hold on

        %% plot pw density CI
        figure(f3)
        errorbar(ax3,experiment_num,[MeanPwDensity(2)], ...
                [MeanPwDensity(2)-MeanPwDensity(1)], ...
                [MeanPwDensity(3)-MeanPwDensity(2)],'o')
        hold on


        %% plot column spacing CI
        figure(f4)
        errorbar(ax4,experiment_num,[MeanSpacing_mm(2)], ...
                [MeanSpacing_mm(2)-MeanSpacing_mm(1)], ...
                [MeanSpacing_mm(3)-MeanSpacing_mm(2)],'o')
        hold on

        
    end

    %% Load and plot data Mouse lemur paper
    file_path = '/home/michael/Cloud/Cloud/PhD/Writing/phd_thesis/Dunnart_Paper/DataMouseLemurPaper/';
    addpath(file_path)
    DataMouseLemurPaper = read_data_tables(file_path);

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.ferret,colorlist(2));

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.cat,colorlist(3));

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.macaque,colorlist(4));

    [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,DataMouseLemurPaper.microcebus,colorlist(5));


    %% Finish Plots and save figures
    LineWidth=3;
    %% save figure1
    figure(f1)
    WeightRange = [1 11*10^3];%[(min(weight_list)*0.9) (max(weight_list)*1.1)];
    plot(ax1,WeightRange,(WeightRange-WeightRange+1)*pi,'--','DisplayName','\pi')

    title(PwDensityType)
    set(gca,'xscale','log')
    set(gca,'DefaultLineLineWidth',LineWidth)
    xlim([1 11*10^3])
    xticks(10.^(0:4))
    ylim([0 6])
    ylabel('Pinwheel Density')
    xlabel('Bodey Weight [g]')
    legend('Location','northeastoutside')


    savefig(f1,[FigureFile '.fig'])
    print(f1,'-depsc2', [FigureFile,'.eps']);

    %% save figure2
    figure(f2)

    %title('Column Spacing [mm]')
    set(gca,'DefaultLineLineWidth',LineWidth)
    set(gca,'xscale','log')
    xlim([1 11*10^3])
    xticks(10.^(0:4))
    ylim([0 1.5])
    ylabel('Column Spacing [mm]')
    xlabel('Bodey Weight [g]')
    legend('Location','northeastoutside')


    savefig(f2,[FigureFileSpacing '.fig'])
    print(f2,'-depsc2', [FigureFileSpacing,'.eps']);

    % %% calc mean corrected pinwheel density BS all animals
    % figure(f3)
    % % MeanPwDensityAnimals = MeanPwDensityAnimals/experiment_Num;
    % %PwDensityCI = bootstrap_ci(SamplesBS,MeanPwDensityAnimals,SamplesJS,alpha);
    % % PwDensityCI(1)=prctile(SamplesBS,2.5);
    % % PwDensityCI(2)=prctile(SamplesBS,100-2.5);
    % % disp(std(SamplesBS)/sqrt(length(SamplesBS))*1.96)
    % BSMeanPwDens=mean(SamplesBS,1); 
    % SortedBSMeanPwDens = sort(BSMeanPwDens);
    % MeanPwDensityAnimals = mean(BSMeanPwDens,'all');
    % PwDensityCI(1)= SortedBSMeanPwDens(floor(length(BSMeanPwDens)*alpha/2));
    % PwDensityCI(2)= SortedBSMeanPwDens(ceil(length(BSMeanPwDens)*(1-alpha/2)));
    % errorbar(ax3,0,[MeanPwDensityAnimals], ...
    %         [PwDensityCI(1)-MeanPwDensityAnimals], ...
    %         [PwDensityCI(2)-MeanPwDensityAnimals],'or')
    % disp([PwDensityCI(1) MeanPwDensityAnimals PwDensityCI(2)])
    % 
    % %% calc mean corrected pinwheel density BS animals
    % MeanPwDensityAnimals = mean(MeanPwDensities,'all');
    % samplesBS = BootstrapSamplesMean1D(MeanPwDensities,N_bootstrapSamples,1);
    % samplesJK = JackknifeSamplesMean1D(MeanPwDensities);
    % PwDensityCI = bootstrap_ci(samplesBS,mean(MeanPwDensities),samplesJK,alpha);
    % errorbar(ax3,0,[MeanPwDensityAnimals], ...
    %         [PwDensityCI(1)-MeanPwDensityAnimals], ...
    %         [PwDensityCI(2)-MeanPwDensityAnimals],'xr')
    % disp([PwDensityCI(1) MeanPwDensityAnimals PwDensityCI(2)])


    % %% calc mean with NumHypercolumns weight
    % weight = NumHypercolumns./mean(NumHypercolumns,'all');
    % 
    % SamplesBS=SamplesBS.*weight;
    % 
    % BSMeanPwDens=mean(SamplesBS,1);
    % SortedBSMeanPwDens = sort(BSMeanPwDens);
    % MeanPwDensityAnimals = mean(BSMeanPwDens,'all');
    % PwDensityCI(1)= SortedBSMeanPwDens(floor(length(BSMeanPwDens)*alpha/2));
    % PwDensityCI(2)= SortedBSMeanPwDens(ceil(length(BSMeanPwDens)*(1-alpha/2)));
    % errorbar(ax3,0,[MeanPwDensityAnimals], ...
    %         [PwDensityCI(1)-MeanPwDensityAnimals], ...
    %         [PwDensityCI(2)-MeanPwDensityAnimals],'+r')
    % disp([PwDensityCI(1) MeanPwDensityAnimals PwDensityCI(2)])

    %% save figure3
    figure(f3)
    hold on
    plot([0 max(experiment_num_list)+1], [pi pi],"DisplayName","\pi")
    title(PwDensityType)
    xlim([0 max(experiment_num_list)+1])
    set(gca,'DefaultLineLineWidth',LineWidth)
    ylim([0 8])
    ylabel('Pinwheel Density')
    xlabel('animal')
    legend('Location','northeastoutside')


    savefig(f3,[FigureFile 'CI.fig'])
    print(f3,'-depsc2', [FigureFile 'CI.eps']);


    %% save figure4
    figure(f4)
    hold on
    title('Column Spacing [mm]')
    xlim([0 max(experiment_num_list)+1])
    set(gca,'DefaultLineLineWidth',LineWidth)
    ylim([0 1.5])
    ylabel('Column Spacing [mm]')
    xlabel('animal')
    legend('Location','northeastoutside')

    savefig(f4,[FigureFileSpacing 'CI.fig'])
    print(f4,'-depsc2', [FigureFileSpacing 'CI.eps']);


    %% plot corr Number of Hypercolumns vs MeanPwDensity
    figure
    plot(NumHypercolumns,MeanPwDensities,'*')
    xlabel('Number Hypercolumns')
    ylabel('mean pw density')
    print('-depsc2', [FigureFile 'PwDensVsHypecolumnN.eps']);

    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

function [f1,ax1,f2,ax2]=PlotDataMouseLemurPaper(f1,ax1,f2,ax2,data,color)

    data.PwD

    figure(f1)
    plot(ax1,data.Weight_g,data.PwD,'o','DisplayName',data.legend_text,'Color',color)
    hold on

    figure(f2)
    plot(ax2,data.Weight_g,data.WL_mm,'o','DisplayName',data.legend_text,'Color',color)
    hold on

end