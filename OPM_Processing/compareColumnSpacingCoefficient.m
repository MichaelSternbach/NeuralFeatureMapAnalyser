function compareColumnSpacingCoefficient(animal,SpecNumList,MainResultFolder,AnimalDataFolder,smallest_w_mm,w_step_mm,largest_w_mm,FilterMap)

    %% set default values
    if nargin<4
        AnimalDataFolder = '~/CIDBN/';
    end
    if nargin <5
        smallest_w_mm = 0.1;
    end
    if nargin < 6
        w_step_mm = 0.1;
    end
    if nargin < 7
        largest_w_mm = 2;
    end
    if nargin < 8
        FilterMap = false;
    end
    getCI = false;

    %% set Figure file name
    mkdir(AnimalDataFolder)
    FigureFileName = [MainResultFolder 'compareColumnSpacingCoefficients_' animal];
    rm_cmd = ['rm -f ' FigureFileName '.ps'];
    disp(rm_cmd)
    system(rm_cmd)

    %% set up figure
    f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
    t = tiledlayout(floor(length(SpecNumList)/3),3, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout

    for specimen_num = SpecNumList


        %% data folder
        ResultDataFolder = [ MainResultFolder lower(animal) '/' lower(animal) num2str(specimen_num) '/'];
        mkdir(ResultDataFolder)

        %% plot parameter
        FigureFileName = [ResultDataFolder 'ResultOPM_' animal num2str(specimen_num) '.eps'];
        border_ColumnSpacing = 20;
        position_ColumnSpacing_txt = [10 10];
        x_hypercolumns = 1;

        %% load data
        disp('load and prepare data...')
        [data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,specimen_num,AnimalDataFolder);
        data_obj.activateGIF(true,4)
        disp('-----------------------')

        %% get column spacing
        [average_spacing_mm,~,~,WavletCoefficient] =  getColumnsSpacing(data_obj,ResultDataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap);
         
        %% plot column spacing coefficients and mean column spacing
        ax = nexttile(t);
        plot(ax,WavletCoefficient.X, WavletCoefficient.Y_mean);
        plot(ax,WavletCoefficient.XI, WavletCoefficient.YI_mean);
        plot(ax,[average_spacing_mm*1000 average_spacing_mm*1000],[min(WavletCoefficient.Y_mean) max(WavletCoefficient.Y_mean)],'r--');
        title([num2str(specimen_num) '. ' data_info.ID])
        xlabel('Scale [mu m]');
        ylabel('Wavelet Coeff');
    end
    
    %% save figure
    print(f,FigureFileName,'-depsc')
end
