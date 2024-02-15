function PlotPwDensityMultiAnimal(PwDensityType,animal_list,experiment_Num,AnimalDataFolder,DataFolderMain,FigureFolder)
    close all
    
    %% parameter
    experiment_Num = checkFormatNum(experiment_Num);
    
%     PwDensityType= 'MeanPwDensity';
    
    %% figure folder
    mkdir(FigureFolder)
    
    FigureFile = [FigureFolder PwDensityType];
    FigureFileSpacing = [FigureFolder 'ColumnSpacing'];
    
%     rm_cmd = ['rm -f ' FigureFile '.ps'];
%     disp(rm_cmd)
%     system(rm_cmd)
    

    %% make Figure
    colorlist = ["r","g","b","c","m","y"];
    f1 = figure();
    f1.Position = [1 1 300 200];
    %set(gca,'LineWidth',4)
    %set(gca,'FontSize',16)
    
    ax1= axes('Parent',f1);
    
    
    
    f2 = figure();
    f2.Position = [1 1 300 200];
    ax2= axes('Parent',f2);
    weight_list = [];
    
    
    
    for ii = 1:length(animal_list)
        animal_ = animal_list(ii);
        animal=convertStringsToChars(animal_);
        if ~isa(experiment_Num,'struct') && (length(experiment_Num) == length(animal_list))
                experiment_Num = 1: experiment_Num(ii);
        else
            experiment_Num_list = experiment_Num.(animal_);
        end
        for experiment_num = experiment_Num_list

            %% data folder
            DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];

            %% load animal info
            [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,1,false,AnimalDataFolder);
            %[data_info,data_path] = info_handle(animal,experiment_num,AnimalDataFolder);

            %% load spacing data
            %[average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,DataFolder,false);
            SpacingFile = [DataFolder 'MapSpacing_' data_info.ID '.mat'];
            load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI')
    %         z = data_obj.filter_map(data_obj.read_map());


            %% load pinwheel data
            PwInfoFile = [DataFolder 'PwInfo_' data_info.ID '.mat'];
            load(PwInfoFile,'PwInfo')

    %         getCI = false;
    %         do_plotting=0;
    %         llp_cutoffs = linspace(0.01, 1,100);
    %         beta=0.5;
    % 
    %         Bootstrapsamples = 100;
    %         data_obj.prepare_samples_array(Bootstrapsamples)
    %         PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta);
    %         
            %% load CI spacing data
            CISpacingFile = [DataFolder 'CI_MapSpacing_' data_info.ID '.mat'];
            load(CISpacingFile,'jackstat_average_spacing_mm','average_spacings_mm')


            %% load pinwheel CI data
            CIPwFile = [DataFolder 'CI_PwDensity_' data_info.ID '.mat'];
            load(CIPwFile,'alpha','PwInfosBS','PwInfosJS')

            %% fix naming for plotting
            data_info.ID = replace(data_info.ID,'_',' ');

            %% get pinwheel density CI     

            MeanPwDensityCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,PwDensityType,false,alpha);
            MeanPwDensity = num2str([MeanPwDensityCI(1) PwInfo.(PwDensityType) MeanPwDensityCI(2)]);
            disp([data_info.ID '  ' MeanPwDensity])

            figure(f1)
            errorbar(ax1,[data_obj.info.weight_in_grams],[PwInfo.(PwDensityType)],[PwInfo.(PwDensityType)-MeanPwDensityCI(1)],[MeanPwDensityCI(2)-PwInfo.(PwDensityType)],'o','DisplayName',data_info.ID,'Color',colorlist(ii))
            hold on

            %% plot spacing CI
            CI_average_spacing_mm = bootstrap_ci(average_spacings_mm(2:end),average_spacing_mm,jackstat_average_spacing_mm,alpha);
            MeanSpacing_mm = {num2str([CI_average_spacing_mm(1) average_spacing_mm CI_average_spacing_mm(2)])};
            disp([MeanSpacing_mm ' mm'])

            figure(f2)
            errorbar(ax2,[data_obj.info.weight_in_grams],[average_spacing_mm],[average_spacing_mm-CI_average_spacing_mm(1)],[CI_average_spacing_mm(2)-average_spacing_mm],'o','DisplayName',data_info.ID,'Color',colorlist(ii))
            hold on

            weight_list =[weight_list data_obj.info.weight_in_grams];
        end
    end
    
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

    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     return 
end







