function PlotResultFactSheets(experiment_num_list,animal,AnimalDataFolder,DataFolderMain,GIF,FigureFolder,BiasDataFolder,width_scale_pix,color_contur,linewidth,Fontsize,alpha)
    
    %PlotResultFactSheets(1:9,'dunnart','~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/OrientationPrefernceMapProcessing/AllData&Results/DataHPC_GIF_adaptedFilter/',4)

    close all

    %% set default values
    if nargin < 5
        GIF = 4;
    end
    if nargin < 6
        FigureFile = [DataFolderMain animal '/' animal 'FactSheet'];
    else
        FigureFile = [FigureFolder animal 'FactSheet'];
    end
    if nargin < 7
        BiasDataFolder = '';
        %%BiasDataFolder = ['~/Cloud/Cloud/PhD/Writing/phd_thesis/OPM_Methods/MakeNoiseFromDataROI_ColumnSpacing/ResultDatav0/' animal '/'];
    end
    if nargin < 8
        width_scale_pix = 15;
    end
    if nargin < 9
        color_contur = [0.6 0.6 0.6];%%'grey';
    end
    if nargin < 10
        linewidth = 0.3;
    end
    if nargin < 11
        Fontsize = 20;
    end

    if nargin < 12
        alpha = 0.3180;%0.05;
    end
    
    if alpha ~= 0.3180
        CI_name = [num2str(round((1-alpha)*100,1)) ' CI '];
    else
        CI_name = 'SE ';
    end

    %% check parameter format
    experiment_num_list = checkFormatNum(experiment_num_list);
    

    
    rm_cmd = ['rm -f ' FigureFile '.ps'];
    disp(rm_cmd)
    system(rm_cmd)

    %% variable to collect CI data
    CI_angle = [];
    CI_Abs = [];
    CI_Pw = [];
    CI_CS = [];
    

    %% loop over experiments
    first = true;
    for experiment_num = experiment_num_list
        
        %% data folder
        DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
        
        %% load animal data
        [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder);
        if GIF>0
            data_obj.activateGIF(true,GIF);
        end
        
        %% plot page
        [AnimalTable,ResultData] = PlotFactSheetPage(data_info,data_obj,DataFolder,FigureFile,...
                                        BiasDataFolder,width_scale_pix,color_contur,linewidth,Fontsize,alpha);
        
        if first
            AllAnimalTable = AnimalTable;
            first = false;
        else
            AllAnimalTable = [AllAnimalTable; AnimalTable];
        end
        
        %% save CI data
        if ~isempty(ResultData)
            CI_angle = [CI_angle; ResultData.CI_angle];
            CI_Abs = [CI_Abs; ResultData.CI_Abs];
            CI_Pw = [CI_Pw ResultData.PwCI];
            CI_CS = [CI_CS; ResultData.CI_CS];
        end

        

    end

    %% plot all data CI CUMPDF

    f = figure();
    t = tiledlayout(2,2);
    title(t,'Data CI CUMPDF All Animals')

    %% plot CPDF CI Maps Angle
    ax = nexttile;
    plotCPDF(CI_angle,'','-',ax)
    title([CI_name 'Prefernce CPDF'])
    xlabel([CI_name 'Prefernce ≤ X [°]'])
    ylabel('% of pix.')
    %xlim([0,1])
    axis(ax,'square')

    %% plot CPDF CI Maps Abs
    ax = nexttile;
    plotCPDF(CI_Abs,'','-',ax)
    title(['rel. ' CI_name 'Selectivity CPDF'])
    xlabel([CI_name 'Selectivity ≤ X [Selectivity]'])
    ylabel('% of pix.')
    %xlim([0,1])
    axis(ax,'square')

    %% plot CPDF CI Maps Pw
    ax = nexttile;
    plotCPDF(CI_Pw,'','-',ax)
    title(['Pinwheel ' CI_name 'Size CPDF'])
    xlabel(['PW ' CI_name 'eff. r ≤ x [mm]'])
    ylabel('% of pinwheels')
    %xlim([0,1])
    axis(ax,'square')

    %% plot CPDF CI Maps CS
    ax = nexttile;
    plotCPDF(CI_CS,'','-',ax)
    title(['rel. ' CI_name 'ColumnSpacing CPDF'])
    xlabel(['rel. ' CI_name '\Lambda/ ≤ X '])
    ylabel('% of pix')
    %xlim([0,1])
    axis(ax,'square')

    %% save figure
    set(gcf,'PaperPositionMode','auto')
    print(f, '-dpsc','-fillpage', '-append', [FigureFile '.ps'])
    
    %% save animal table to latex and seperate page
%     AllAnimalTable = rows2vars(AllAnimalTable);
%     AllAnimalTable.Properties.VariableNames(1) = "Name";
    disp(AllAnimalTable)
    table2latex(AllAnimalTable, [FigureFile 'Table.tex'])
%     FigureAllAnimalTable = plotTable(AllAnimalTable);
%     print(FigureAllAnimalTable, '-dpsc','-fillpage', '-append', [FigureFile '.ps'])
    
    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%ps2pdf dunnart_FactSheet.ps dunnart_FactSheet.pdf
    system(['ps2pdf ' FigureFile '.ps ' FigureFile '.pdf'])
    return 
end




