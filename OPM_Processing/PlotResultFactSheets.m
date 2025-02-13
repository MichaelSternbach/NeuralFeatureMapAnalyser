function PlotResultFactSheets(experiment_num_list,animal,AnimalDataFolder,DataFolderMain,GIF,FigureFolder,BiasDataFolder,width_scale_pix,color_contur,linewidth,Fontsize)
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
    
    %% check parameter format
    experiment_num_list = checkFormatNum(experiment_num_list);
    

    
    rm_cmd = ['rm -f ' FigureFile '.ps'];
    disp(rm_cmd)
    system(rm_cmd)
    

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
        AnimalTable = PlotFactSheetPage(animal,experiment_num,data_info,data_obj,DataFolder,FigureFile,...
                                        BiasDataFolder,width_scale_pix,color_contur,linewidth,Fontsize);
        
        if first
            AllAnimalTable = AnimalTable;
            first = false;
        else
            AllAnimalTable = [AllAnimalTable; AnimalTable];
        end
        

        

    end
    
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





