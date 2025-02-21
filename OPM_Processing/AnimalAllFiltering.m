function AnimalAllFiltering(animal,experiment_num_list,AnimalDataFolder,DataFolderMain,...
    GIF_SNTH,FigureFolder,power_spectrum_range_mm,lowpass_cutoffs_mm,MinLengthPlateau_mm,sigma)
    
    close all
    %% parameter
    experiment_num_list = checkFormatNum(experiment_num_list);
    if nargin <5
        GIF_SNTH = 4;
    end
    if nargin<6
        FigureFolder = DataFolderMain;
    end
    if nargin<7
        power_spectrum_range_mm = 0.2:0.01:2.5;
    end
    if nargin<8
        lowpass_cutoffs_mm =0.1:0.01:0.4;
    end
    if nargin <9
        MinLengthPlateau_mm = 0.1;
    end
    if nargin <10
        sigma=0.1;
    end
    
    
    %% figure folder and files
    if ~isfolder(FigureFolder)
        mkdir(FigureFolder) 
    end

    FigureFile1 = [FigureFolder animal '_PowerSpectrums'];
    FigureFile2 = [FigureFolder animal '_PowerProfiles'];
    FigureFile3 = [FigureFolder animal '_PwPlateaus'];
    FigureFile4 = [FigureFolder animal '_PwDensPlateaus'];
    
    rm_cmd = ['rm -f ' FigureFile1 '.ps'];
    system(rm_cmd)
    rm_cmd = ['rm -f ' FigureFile2 '.ps'];
    system(rm_cmd)
    rm_cmd = ['rm -f ' FigureFile3 '.ps'];
    system(rm_cmd)
    rm_cmd = ['rm -f ' FigureFile4 '.ps'];
    system(rm_cmd)
    
    
    
    %% make Figure
    N_Animals = length(experiment_num_list);

    f1 = figure(1);
    t1 = tiledlayout(ceil(N_Animals/3),3);
    title(t1,'Power Spectra')
    
    f2 = figure(2);
    t2 = tiledlayout(ceil(N_Animals/3),3);
    title(t2,'Power Profiles')

    f3 = figure(3);
    t3 = tiledlayout(ceil(N_Animals/3),3);
    title(t2,'Pinwheel Number Plateaus')
    
    f4 = figure(4);
    t4 = tiledlayout(ceil(N_Animals/3),3);
    title(t2,'Pinwheel Density Plateaus')

    for experiment_num = experiment_num_list
        
        %% data folder
        DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
        
        %% animal data
        [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder);
        if GIF_SNTH>0
            data_obj.activateGIF(true,GIF_SNTH)
        end
        z = data_obj.read_map();
        
        %% column spacing
%         [average_spacing_mm,local_spacing_mm,~] =  getColumnsSpacing(data_obj,DataFolder,false);
        SpacingFile = [DataFolder 'MapSpacingFiltered_' data_obj.info.ID '.mat'];
        load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI','WavletCoefficient')

        
        %% get power profile data
        power_profile = define_filter_settings(data_info,data_obj.ROI,z,power_spectrum_range_mm);

        %% get Pw plateu data
        PwResultData = findLowpassPwNumbers(data_obj,lowpass_cutoffs_mm,true,local_spacing_mm,sigma);
        

        %% plot powerspectrum OPM
        figure(1);
        nexttile(t1);
        plot(1./power_profile.scale_mm,power_profile.values_kspace);
        hold on
        plot([1./data_obj.filter_parameters.lowpass 1./data_obj.filter_parameters.lowpass],[min(power_profile.values_kspace,[],'all') max(power_profile.values_kspace,[],'all')],'DisplayName','Lowpass Cutoff')
        hold on
        plot([1./data_obj.filter_parameters.highpass 1./data_obj.filter_parameters.highpass],[min(power_profile.values_kspace,[],'all') max(power_profile.values_kspace,[],'all')],'DisplayName','Highpass Cutoff')
        xlabel('spatial frequwncy [1/mm]')
        ylabel('Power')
        axis square
        title(data_obj.info.ID)

        %% plot powerspectrum OPM
        figure(2);
        nexttile(t2);
        plot(power_profile.scale_mm,power_profile.values);
        hold on
        plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
        hold on
        plot([data_obj.filter_parameters.highpass data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
        xlabel('spatial lengthscale [mm]')
        ylabel('Power')
        axis square
        title(data_obj.info.ID)

        %% plot Pw density plateau
        figure(3);
        nexttile(t3);
        plotPlateau(PwResultData.lowpass_cutoffs,PwResultData.NumberPw,data_obj.filter_parameters.lowpass,MinLengthPlateau_mm)
        title(data_obj.info.ID)
        legend()
        xlabel('lowpass cutoff [mm]')
        ylabel('Pinwheel Number')
        axis square;    
        
        %% plot Pw density plateau
        figure(4);
        nexttile(t4);
        plotPlateau(PwResultData.lowpass_cutoffs,PwResultData.DensityPw,data_obj.filter_parameters.lowpass,MinLengthPlateau_mm)
        title(data_obj.info.ID)
        legend()
        xlabel('lowpass cutoff [mm]')
        ylabel('Pinwheel Density')
        axis square;    
        
    end
    

    print(f1, '-dpsc','-fillpage', [FigureFile1 '.ps'])%'-append',
    print(f2, '-dpsc','-fillpage', [FigureFile2 '.ps'])%'-append',
    print(f3, '-dpsc','-fillpage', [FigureFile3 '.ps'])%'-append',
    print(f4, '-dpsc','-fillpage', [FigureFile4 '.ps'])%'-append',
    
    
    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return 
end

function plotPlateau(x,y,MarkerX,MinLengthPlateau,FitRange)

    %% input
    if nargin <4
        MinLengthPlateau = (FitRange(2)-FitRange(1))*0.1;
    end
    if nargin <5
        FitRange = [min(x) max(x)];
    end
    
    %% fit plateu
    
    
    %% plot data
    plot(x,y,'x','DisplayName','data')
    hold on

    %% plot lowpass cutoff
    minY = min(y)*0.9;
    maxY = max(y)*1.1;
    if isnumeric(MarkerX)
        plot([MarkerX MarkerX],[minY maxY],'b')
    end

    %% add plateau fit
    try
        [~, plt, param, ssel, ~]= fit_piecewise_linear_Manuel(x,y,FitRange,MinLengthPlateau,false);
        fitY = PlateauFit(x,plt, param, ssel);
        hold on;
        plot(x,fitY,'r','DisplayName','fit')
    catch
        disp('Plateau could not be found!')
    end


    %% add ylim
    xlim(FitRange)
    ylim([minY maxY])

end

function yy = PlateauFit(x,plt, params, ssel)
    %plt = [mnplt_min mxplt_min]
    mnplt_min = plt(1);
    mxplt_min = plt(2);
    % Compute the actual fit function
    if isequal(ssel, [1, 1, 1])
        yy = params(1) + (-x+mnplt_min).*((-x+mnplt_min)>0)*params(2) + (-x+mxplt_min).*((-x+mxplt_min)<0)*params(3);       
    elseif isequal(ssel, [1, 0, 1])
        yy = params(1) +  (-x+mxplt_min).*((-x+mxplt_min)<0)*params(2);
    elseif isequal(ssel, [1, 1, 0])
        yy = params(1) + (-x+mnplt_min).*((-x+mnplt_min)>0)*params(2);
    end

end



