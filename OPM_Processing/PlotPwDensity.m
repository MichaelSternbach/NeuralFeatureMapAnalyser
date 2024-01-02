function PlotPwDensity(PwDensityType,animal,experiment_Num,AnimalDataFolder,DataFolderMain,FigureFolder)
    close all
    
    %% parameter
    experiment_Num = checkFormatNum(experiment_Num);
    
%     PwDensityType= 'MeanPwDensity';
    
    %% figure folder
    mkdir(FigureFolder)
    
    FigureFile = [FigureFolder animal '_' PwDensityType];
    FigureFileSpacing = [FigureFolder animal '_ColumnSpacing'];
    
    rm_cmd = ['rm -f ' FigureFile '.ps'];
    disp(rm_cmd)
    system(rm_cmd)
    

    %% make Figure

    f1 = figure();
    ax1= axes('Parent',f1);
    
    f2 = figure();
    ax2= axes('Parent',f2);

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
        
        %% load CI spacing data
        CISpacingFile = [DataFolder 'CI_MapSpacing_' data_obj.info.ID '.mat'];
        load(CISpacingFile,'CI_average_spacing_mm')
        

        %% load pinwheel CI data
        CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
        load(CIPwFile,'alpha','PwInfosBS','PwInfosJS')

        %% get pinwheel density CI     
        
        MeanPwDensityCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,PwDensityType,false,alpha);
        MeanPwDensity = num2str([MeanPwDensityCI(1) PwInfo.(PwDensityType) MeanPwDensityCI(2)]);
        disp([data_info.ID '  ' MeanPwDensity])
        
        figure(f1)
        errorbar(ax1,[experiment_num],[PwInfo.(PwDensityType)],[PwInfo.(PwDensityType)-MeanPwDensityCI(1)],[MeanPwDensityCI(2)-PwInfo.(PwDensityType)],'o','DisplayName',data_info.ID)
        hold on
   
        %% plot spacing CI
        MeanSpacing_mm = {num2str([CI_average_spacing_mm(1) average_spacing_mm CI_average_spacing_mm(2)])};
        disp([MeanSpacing_mm ' mm'])
        
        figure(f2)
        errorbar(ax2,[experiment_num],[average_spacing_mm],[average_spacing_mm-CI_average_spacing_mm(1)],[CI_average_spacing_mm(2)-average_spacing_mm],'o','DisplayName',data_info.ID)
        hold on

        
    end
    
    %% save figure1
    figure(f1)
    animalRange = 0:experiment_Num+1;
    plot(ax1,animalRange,(animalRange-animalRange+1)*pi,'DisplayName','\pi')
    
    title(PwDensityType)
    set(gca,'DefaultLineLineWidth',2)
    xlim([0 experiment_Num+1])
    ylabel('Pinwheel Density')
    xlabel('Animal')
    legend('Location','northeastoutside')
    
    
    savefig(f1,[FigureFile '.fig'])
    print(f1, '-dpsc','-fillpage', '-append', [FigureFile '.ps'])
    
    %% save figure2
    figure(f2)

    %title('Column Spacing [mm]')
    set(gca,'DefaultLineLineWidth',2)
    xlim([0 experiment_Num+1])
    ylabel('Column Spacing [mm]')
    xlabel('Animal')
    legend('Location','northeastoutside')
    
    
    savefig(f2,[FigureFileSpacing '.fig'])
    print(f2, '-dpsc','-fillpage', '-append', [FigureFileSpacing '.ps'])

    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     return 
end







