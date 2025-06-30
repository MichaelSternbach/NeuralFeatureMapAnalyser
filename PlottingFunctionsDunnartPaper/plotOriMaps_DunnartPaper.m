function plotOriMaps_DunnartPaper(experiments,AnimalDataFolder,DataFolderMain,...
    GIF_SN_TH,data_info_file_list,smallest_w_mm,largest_w_mm,w_step_mm,Bootstrapsamples)
    
%     experiments.dunnart = [1 2 4 5 6 7 8 9];
%     addpath /home/michael/Cloud/Cloud/PhD/MarsupialData/OrientationPrefernceMapProcessing/OPM_Processing
%     plotOriMaps_DunnartPaper(experiments,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/OrientationPrefernceMapProcessing/AllData&Results/DataHPC_GIF_adaptedFilter/',4,'/home/michael/Cloud/Cloud/PhD/MarsupialData/OrientationPrefernceMapProcessing/OPM_Processing/experiment_info.csv',0.1,2,0.05,100)


    close all
    %% parameter
    %experiment_num_list = checkFormatNum(experiment_num_list);
    if nargin <5
        GIF_SN_TH = 4;
    end
    if nargin < 6 %|| isempty(data_info_file_list)
        function_dir = fileparts(mfilename('fullpath')); % Get function directory
        data_info_file_list = fullfile(function_dir, 'experiment_info.csv');
    end

    rotate = -90;
    linewidth = 1;
    %Bootstrapsamples = 100;
    color_contur = [0.6 0.6 0.6];%%'grey';
%     
%     %% set filter settings
%     lowpass_mm = 0.25;
%     highpass_mm = 0.56;
    
    %% figure folder and files
    if ~isfolder(DataFolderMain)
        mkdir(DataFolderMain) 
    end

    %% set alpha
%     alpha = 0.05;
    alpha = 0.3180; % SE

    %% make data folder
    if ~isfolder(DataFolderMain)
        mkdir(DataFolderMain) 
    end
  
    %% extract experiment_list from keys of the struct
    experiment_list = fieldnames(experiments);
    result_data = struct();

    for ii = 1:length(experiment_list)
        experiment = experiment_list{ii};
        experiment_num_list = experiments.(experiment);
        result_data.(experiment) = cell(1,max(experiment_num_list));
        DataFolderExperiment = [DataFolderMain experiment '/'];
        if ~isfolder(DataFolderExperiment)
            mkdir(DataFolderExperiment) 
        end
        for idx = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(idx);
            %% make datafolder
            DataFolder = [DataFolderExperiment experiment num2str(experiment_num) '/'];
            if ~isfolder(DataFolder)
                mkdir(DataFolder) 
            end

            %% plot one animal and extract data
            [result_data.(experiment){experiment_num}] = plotOneAnimal(experiment,experiment_num,AnimalDataFolder,DataFolder,data_info_file_list,GIF_SN_TH,smallest_w_mm,largest_w_mm,w_step_mm,Bootstrapsamples,alpha);
        end

        %% comparison plots one species

        %% SE distributions CPDF
        all_CI_dir = [];
        f1 = figure;
        f1.Position = [1 1 200 200];
        for jj = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(jj);
            data = result_data.(experiment){experiment_num};
            CI_dir = data.CI.BCA.CI_angle;
            ROI = data.ROI;
            CI_dir = CI_dir(ROI);
            all_CI_dir = [all_CI_dir; CI_dir];
            plot(sort(CI_dir),linspace(0,1,length(CI_dir)),'DisplayName',[experiment num2str(experiment_num)])
            hold on
        end
        plot(sort(all_CI_dir),linspace(0,1,length(all_CI_dir)),'--black',linewidth=2)
        xlabel('SE [°]')
        ylabel('Cumulative Probability')
        print(f1,'-dpsc', [DataFolderMain experiment 'SE_CDF.ps'])
        %print(f1, [DataFolderMain experiment 'SE_CDF'], '-depsc', '-vectors');

        %% SizesCI distribution CPDF
        all_PwCI = [];
        f1 = figure;
        f1.Position = [1 1 200 200];
        for jj = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(jj);
            data = result_data.(experiment){experiment_num};
            PwCI = data.PwCI_r*1000;
            all_PwCI = [all_PwCI PwCI];
            plot(sort(PwCI),linspace(0,1,length(PwCI)),'DisplayName',[experiment num2str(experiment_num)])
            hold on
        end
        plot(sort(all_PwCI),linspace(0,1,length(all_PwCI)),'--black',linewidth=2)
        xlabel('Pinwheel SE eff. r [µm]')
        ylabel('Cumulative Probability')
        print(f1,'-dpsc', [DataFolderMain experiment 'SizesCI_CDF.ps'])


        %% violine plot peaks_test comparioson
        peaks_test_list = cell(1,length(experiment_num_list));
        peaks_test_rand_list = cell(1,length(experiment_num_list));
        x_labels = cell(1,length(experiment_num_list));
        for jj = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(jj);
            data = result_data.(experiment){experiment_num};
            peaks_test_list{jj} = data.peaks_test/max(data.peaks_test);
            peaks_test_rand_list{jj} = data.peaks_test_rand/max(data.peaks_test);
            x_labels{jj} = experiment_num;
        end

        f1 = figure;
        f1.Position = [1 1 70*length(peaks_test_list) 250];
        violin(peaks_test_list,'facecolor',[0 0 1],'xlabel',x_labels);
        hold on;
        violin(peaks_test_rand_list,'facecolor',[1 0 0],'xlabel',x_labels);
        ylabel('modularity index')
        ylim([0 1.5])
        xlabel(experiment)
        print(f1,'-dpsc', [DataFolderMain experiment 'modularity_peaks.ps'])

        %% violine plot pinwheel stability indicies comparioson
        stability_indicies_list = cell(1,length(experiment_num_list));
        stability_indicies_rand_list = cell(1,length(experiment_num_list));
        x_labels = cell(1,length(experiment_num_list));
        for jj = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(jj);
            data = result_data.(experiment){experiment_num};
            stability_indicies_list{jj} = data.StabilityIndicies.data;
            stability_indicies_rand_list{jj} = data.StabilityIndicies.data_rand;
            x_labels{jj} = experiment_num;
        end
        f1 = figure;
        f1.Position = [1 1 70*length(stability_indicies_list) 250];
        violin(stability_indicies_list,'facecolor',[0 0 1],'xlabel',x_labels);
        hold on;
        violin(stability_indicies_rand_list,'facecolor',[1 0 0],'xlabel',x_labels);
        ylabel('stability index')
        title('test pinwheel configuration')
        ylim([0 1.2])
        xlabel(experiment)
        print(f1,'-dpsc', [DataFolderMain experiment 'pinwheel_stability.ps'])


        %% plot mean spacing vs std local spacing
        mean_spacing_list = [];
        std_local_spacing_list = [];    
        for jj = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(jj);
            data = result_data.(experiment){experiment_num};
            mean_spacing_list = [mean_spacing_list data.average_spacing_mm];
            std_local_spacing_list = [std_local_spacing_list std(data.local_spacing_mm(data.ROI))];
        end

        f1 = figure;
        f1.Position = [1 1 200 200];
        plot(mean_spacing_list,std_local_spacing_list,'o')
        xlabel('mean spacing [mm]')
        ylabel('std local spacing [mm]')
        print(f1,'-dpsc', [DataFolderMain experiment 'mean_std_spacing.ps'])

    
    end


    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return 
end

function [result_data] = plotOneAnimal(animal,experiment_num,AnimalDataFolder,DataFolder,data_info_file_list,GIF_SN_TH,smallest_w_mm,largest_w_mm,w_step_mm,Bootstrapsamples,alpha)
    
    rotate = -90; 
    
    %% get data
    [data_info,data_path,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder,data_info_file_list);
    if GIF_SN_TH>0
        data_obj.activateGIF(true,GIF_SN_TH)
    end
    data_obj.prepare_samples_array(Bootstrapsamples)

    if animal == "dunnart"
        %% set filter settings
        lowpass_mm = 0.25;
        highpass_mm = 0.56;

        data_info.settings.lowpass_mm = lowpass_mm;
        data_info.settings.highpass_mm = highpass_mm;

        data_obj.set_filter_parameters("lowpass",lowpass_mm);
        data_obj.set_filter_parameters("highpass",highpass_mm);
    end


    %% get OPM map
    result_data.z = data_obj.filter_map(data_obj.read_map(),false);

    %% save ROI
    result_data.ROI = data_obj.ROI;
    %% get CI column spacing Ori
    %[average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap,alpha)   
    [result_data.average_spacing_mm,result_data.local_spacing_mm,newROI,~,result_data.CI_average_spacing_mm,result_data.CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,true,true,alpha,false);

    %% get pw infos
    do_plotting=0;
    PwDensitCalcSteps_lowpass_cutoffs_mm = 0.1:0.01:0.6;
    SizeGaussKernelPwDensityCalc = 0.5;
    result_data.PwInfo = getPinwheelInfos(data_obj,result_data.local_spacing_mm,DataFolder,newROI,true,do_plotting,PwDensitCalcSteps_lowpass_cutoffs_mm,SizeGaussKernelPwDensityCalc);
    
    %% get CI
    DoFilter = true;
    result_data.CI = calcCIs(data_obj,alpha,DoFilter,DataFolder);
    %CI.BCA.CI_angle,CI.BCA.CI_Abs,CI.BCA.ROI


    %% make map borders ROI
    rel_x_length = 2;
    rel_y_length = 3;
    [x_range,y_range] = makeRectangleROI(data_obj,result_data.average_spacing_mm,rel_x_length,rel_y_length);

    %% set electrode position
    %compareParameterGIF('dunnart',7,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/',[92 110])
    %close all; compareParameterGIF('dunnart',6,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/',[67 96])
    if animal == "dunnart"
        if experiment_num == 6
            ElectrodePosition = [67 96];
        elseif experiment_num == 7
            ElectrodePosition = [92 110];
        else
            ElectrodePosition = [];
        end
    else
        ElectrodePosition = [];
    end


    %% set all variables
    z = result_data.z;
    ROI = data_obj.ROI;
    PwInfo = result_data.PwInfo;
    CI_angle = result_data.CI.BCA.CI_angle;
    average_spacing_mm=result_data.average_spacing_mm;
%     ROI(ROI)=1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% plot polar map  + scale bar 0.5 mm + CI Contour 
    f1 = figure();
    set(f1, 'InvertHardcopy', 'off'); 
    plot_map(z,ones(size(z)),3*sqrt(mean(abs(z(ROI)).^2)),1)
    title('orientation pref. polar map ')

    hold on
    Xmin = x_range(1);
    Ymin = y_range(1);
    mm = 0.2;
    width_scale_pix = 5;
    spacing_pix = mm * data_info.pix_per_mm;
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-white','LineWidth',2.5)


    hold on; [C, h] = contour(CI_angle,[20 20],'white','LineWidth',2,'LineStyle',':');
    %clabel(C,h,'FontSize',10,'Color','white');

    hold on; [C, h] = contour(CI_angle,[80 80],'white','LineWidth',2,'LineStyle','-');
    %clabel(C,h,'FontSize',10,'Color','white');

    %hold on; contour(ROI,[1 1],'white','linewidth',2.5,'LineStyle','-')
    
    % Create mask for everything outside ROI
    outside_roi = ~ROI;
    % Overlay white semi-transparent shade
    hold on
    h = imagesc(cat(3, ones(size(outside_roi)), ones(size(outside_roi)), ones(size(outside_roi))));  % white RGB image
    set(h, 'AlphaData', 0.7 * outside_roi);  % 0.7 transparency outside ROI, 0 inside

    if ~isempty(ElectrodePosition)
    hold on;
    plot(ElectrodePosition(1),ElectrodePosition(2),'xwhite','MarkerSize',15,'LineWidth',10)
    end
    
    xlim(x_range)
    ylim(y_range)
    hold on; set(gca,'view',[rotate rotate])
    print(f1,'-dpdf', '-image', '-fillpage', '-r1200', ...
     [DataFolder animal num2str(experiment_num) 'PolarOri.pdf'])



    %% plot angle map  + scale bar Lambda + Pw CI 
    f1 = figure();
    plot_map(z,ones(size(z)),0,1)
    xlim(x_range)
    ylim(y_range)
    title('Orientation pref. map ')

    hold on
    Xmin = x_range(1);
    Ymin = y_range(1);
    mm = result_data.average_spacing_mm;
    width_scale_pix = 5;
    spacing_pix = mm * data_info.pix_per_mm;
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-white','LineWidth',2.5)

    hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','MarkerSize',5)
    hold on; SizesCI = getConfidenceRegionPw(PwInfo.pinwheel_stats,data_info.field_size_pix,1-alpha,true,false,false,'white',2);
    result_data.PwCI_r = sqrt(SizesCI/(data_info.pix_per_mm)^2/pi);

%     hold on; contour(CI_angle,[30 30],'white','linewidth',2,'LineStyle','-')

    % Create mask for everything outside ROI
    outside_roi = ~ROI;
    % Overlay white semi-transparent shade
    hold on
    h = imagesc(cat(3, ones(size(outside_roi)), ones(size(outside_roi)), ones(size(outside_roi))));  % white RGB image
    set(h, 'AlphaData', 0.7 * outside_roi);  % 0.7 transparency outside ROI, 0 inside
    

    hold on; set(gca,'view',[rotate rotate])
    print(f1,'-dpdf', '-image', '-fillpage', '-r1200', ...
    [DataFolder animal num2str(experiment_num) 'OriContourCI.pdf'])



    %% plot angle map  + Pw sign 
    f1 = figure();
    %f1.Position = [1 1 100 100];
    plot_map(z,ROI,0,1)
    xlim(x_range)
    ylim(y_range)
    title('pinwheel signs ')
    hold on

    % Plot filled white circles
    scatter(PwInfo.PWxList,PwInfo.PWyList, 100, 'w', 'filled');  % size 100, color white, filled
    hold on;
    


    % Plot white '+' signs over the circles
    idx_pos = PwInfo.signList>0;
    x_p = PwInfo.PWxList(idx_pos);
    y_p = PwInfo.PWyList(idx_pos);
    plot(x_p, y_p, 'black+', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on

    % Plot white '-' signs over the circles
    idx_neg = PwInfo.signList<0;
    x_n = PwInfo.PWxList(idx_neg);
    y_n = PwInfo.PWyList(idx_neg);
    plot(x_n, y_n, 'k_', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on

    
    Xmin = x_range(1);
    Ymin = y_range(1);
    mm = result_data.average_spacing_mm;
    width_scale_pix = 5;
    spacing_pix = mm * data_info.pix_per_mm;
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-white','LineWidth',2.5)


%     hold on; set(gca,'view',[rotate rotate])
    print(f1,'-dpsc', [DataFolder animal num2str(experiment_num) 'OriPwSign.ps'])


    %% plot local_spacing_mm + zero contours + Pws

    f1 = figure();
    
%     plot_mapAbs(result_data.local_spacing_mm,'local spacing [mm]',average_spacing_mm*1.2,average_spacing_mm*.8,data_obj.ROI,gca)
    imagesc(gca,result_data.local_spacing_mm,[round(average_spacing_mm*.7,1) round(average_spacing_mm*1.3,1)]);
    colormap(gca, 'gray');%gray,'jet'
    %colormap hot;
    
    

    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    hold on
    %colorbar('Ticks',[0,0.5,3])
    cb = colorbar(gca) ;
    
    hold on; contour(ROI,[1 1],'white','LineWidth',2)

    hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','MarkerSize',5)

    hold on; contour(real(z.*ROI),[0 0],'LineWidth',2,'Color','white')
    hold on; contour(imag(z.*ROI),[0 0],'LineWidth',2,'Color','white')

    xlim(x_range)
    ylim(y_range)

    

    ax = gca;
    ax_pos = ax.Position;        % [x y width height]
    cb.Position(2) = ax_pos(2);  % align vertical position
    cb.Position(4) = ax_pos(4); 

    hold on; set(gca,'view',[rotate rotate])
    title('local spacing [mm]')
    

    print(f1,'-dpsc', [DataFolder animal num2str(experiment_num) 'local_spacing.ps'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% test modularity
    profile_range_mm = smallest_w_mm:w_step_mm:largest_w_mm;
    [result_data.peaks_test,result_data.peaks_test_rand,result_data.peak_position_mm,result_data.power_profiles,result_data.power_profiles_rand,~,~] = testModularityOPM(data_obj,DataFolder,result_data.average_spacing_mm,profile_range_mm,Bootstrapsamples,false,false);

    %% calc stability indicies
    %result_data.StabilityIndicies=getPwConfigStabilityIndices(data_obj,Bootstrapsamples);
    [~,~,~,norm_sel,norm_sel_rand] = testPWsOPM(data_obj,PwInfo.pinwheel_stats,Bootstrapsamples,DataFolder,false,false);
%     mean_rand = nanmean(norm_sel_rand);
%     result_data.StabilityIndicies.data = 1-abs(norm_sel-mean_rand)/mean_rand;
%     result_data.StabilityIndicies.data_rand = 1-abs(norm_sel_rand-mean_rand)/mean_rand;

    max_ = max([norm_sel_rand norm_sel]);
    result_data.StabilityIndicies.data = 1-(norm_sel)/max_;
    result_data.StabilityIndicies.data_rand = 1-(norm_sel_rand)/max_;

    %[result_data.StabilityIndicies.data,result_data.StabilityIndicies.data_rand] = testPwConfig(data_obj,Bootstrapsamples,DataFolder);
end

function [x_range,y_range] = makeRectangleROI(data_obj,average_spacing_mm,rel_x_length,rel_y_length)

    %% get ROI borders
    [x_range,y_range] = getRangeXY_ROI(data_obj.ROI);

    %% get average spacing pix
    average_spacing_px = average_spacing_mm*data_obj.info.pix_per_mm; 

    %% get new x/y_range at the center of the ROI
    x_center = round(x_range(1) + (x_range(2)-x_range(1))/2);
    y_center = round(y_range(1) + (y_range(2)-y_range(1))/2);
    width_y_range = rel_y_length*average_spacing_px;
    width_x_range = rel_x_length*average_spacing_px;
    
    x_range = [x_center-width_x_range/2 x_center+width_x_range/2];
    y_range = [y_center-width_y_range/2 y_center+width_y_range/2];

end


function [result_data,data_info,data_obj,data_path] = getOriResultData(animal,experiment_num,AnimalDataFolder,data_info_file_list,smallest_w_mm,largest_w_mm,w_step_mm,GIF_SN_TH,Bootstrapsamples,alpha)%smallest_w_mm,largest_w_mm,w_step_mm

    %% get data
    [data_info,data_path,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder,data_info_file_list);
    if GIF_SN_TH>0
        data_obj.activateGIF(true,GIF_SN_TH)
    end
    data_obj.prepare_samples_array(Bootstrapsamples)

    if animal == "dunnart"
        %% set filter settings
        lowpass_mm = 0.25;
        highpass_mm = 0.56;

        data_info.settings.lowpass_mm = lowpass_mm;
        data_info.settings.highpass_mm = highpass_mm;

        data_obj.set_filter_parameters("lowpass",lowpass_mm);
        data_obj.set_filter_parameters("highpass",highpass_mm);
    end


    %% get OPM map
    result_data.z = data_obj.filter_map(data_obj.read_map());

    %% sve ROI
    result_data.ROI = data_obj.ROI;

%     %% get column spacing
%     %use [average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
%     [result_data.average_spacing_mm,result_data.local_spacing_mm,result_data.newROI,result_data.WavletCoefficient] = get_column_spacingManuel(result_data.z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);

    %% CI Ori
    [result_data.CI,~,~] = getCI(data_obj,alpha,'bca',true,false,true);


    % %% get Pw infos
    % lowpass_cutoffs_mm = 0.1:0.01:0.5;
    % SizeGaussKernelPwDensityCalc = 0.5;
    % do_plotting = false;
    % result_data.PwInfo = getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI_PwCS,do_plotting,lowpass_cutoffs_mm,SizeGaussKernelPwDensityCalc);

    %% find pws
    [result_data.PwInfo.NumberPw,result_data.PwInfo.aniso,result_data.PwInfo.x_angle,result_data.PwInfo.PWxList,result_data.PwInfo.PWyList,result_data.PwInfo.signList, result_data.PwInfo.contours] = find_pinwheels(result_data.z,0,data_obj.ROI);
    %result_data.PwInfo.pinwheel_stats = getPinwheelStats(data_obj,result_data.PwInfo.x_angle,result_data.PwInfo.y_angle,result_data.PwInfo.x_list,result_data.PwInfo.y_list);

end

function [result_data,data_obj_dir]=getDirResultData(data_info,data_path,data_obj,smallest_w_mm,largest_w_mm,w_step_mm,GIF_SN_TH,Bootstrapsamples,alpha)
       %% get dir map
       [~,~,DirectionData] = getDirectionData(data_info,data_path,0);

        DirectionData  = real(DirectionData);
    
        %% make data obj dir data
        data_info_dir = data_info;
        stimDir = [data_info.stim_order(find(~isnan(data_info.stim_order))) data_info.stim_order(find(~isnan(data_info.stim_order)))+180];
        stimDir = [stimDir NaN];
        data_info_dir.stim_order=stimDir;
        data_obj_dir = data_handle_corrected(data_info_dir,DirectionData,data_obj.ROI);
        if GIF_SN_TH>0
            data_obj_dir.activateGIF(true,GIF_SN_TH)
        end
        data_obj_dir.prepare_samples_array(Bootstrapsamples)
    
    
        %% get dir map 
        result_data.z_dir = data_obj_dir.filter_map(data_obj_dir.read_map(1,false,true));

        %% save ROI 
        result_data.ROI = data_obj_dir.ROI;
    
        %% calc CI dir and Ori map
        [result_data.CI,~,~] = getCI(data_obj_dir,alpha,'bca',true,true,true);

%         %% calc column spacing
%         %use [average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
%         [result_data.average_spacing_mm,result_data.local_spacing_mm,result_data.newROI,result_data.WavletCoefficient] = get_column_spacingManuel(result_data.z_dir,data_obj_dir.ROI,data_obj_dir.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);

end


%  %% Plot Bloodvessels
%  PlotBloodVessels(BloodVesselImg,ones(size(BloodVesselImg)),1)
%  hold on
%  contour(ROI,[1 1],'red','linewidth',2)
%  hold on
%  mm = 0.2;
%  width_scale_pix = 5;
%  spacing_pix = mm * data_info.pix_per_mm;
%  plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-red')
%  hold on
%  text(Xmin+width_scale_pix,Ymin+width_scale_pix+2,[num2str(mm) ' mm'],'Color','red')
%  set(gca,'xtick',[])
%  set(gca,'ytick',[])
%  xlim(x_range)
%  ylim(y_range)

%  hold on; set(gca,'view',[rotate rotate])
%  print('-depsc', [DataFolder 'BloodVessels_' animal num2str(experiment_num) '.eps'])
%  close

%  %% plot polar map
%  f1 = figure();
%  plot_map(z,data_obj.ROI,3*sqrt(mean(abs(z(ROI)).^2)),1)
%  xlim(x_range)
%  ylim(y_range)
%  title('orientation pref. map ')
%  hold on; set(gca,'view',[rotate rotate])
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'PolarOri.ps'])

%  %% plot ori map big
%  f1=figure;
%  plot_map(z,data_obj_dir.ROI,0,1)

%  hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','MarkerSize',5)
%  hold on; SizesCI = getConfidenceRegionPw(PwInfo.pinwheel_stats,data_info.field_size_pix,1-alpha,true,false,false,'black',2);
%  %getConfidenceRegionPw(pinwheel_stats,field_size_pix,Confidence,plot_CI,Plot_CI_Size_Scale,plotProb,color,linewidth)

%  hold on; contour(NoiseMeasure.CI_ori,[30 30],'white','linewidth',2,'LineStyle','-')
%  %hold on; contour(NoiseMeasure.CI_ori,[40 40],'white','linewidth',2,'LineStyle','--')

%  xlim(x_range)
%  ylim(y_range)
%  title('orientation pref. map ')
%  hold on; set(gca,'view',[rotate rotate])
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'Ori.ps'])

%  %% plot ori map Abs big
%  f1=figure;
%  plot_Abs(abs(z).*ROI,gca)
%  xlim(x_range)
%  ylim(y_range)
%  title('selectivity orientation map')
%  hold on; set(gca,'view',[rotate rotate])
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'OriAbs.ps'])


%  %% plot Pw CI
%  f1 = figure;
%  sigma = 0.05;
%  PWxList = PwInfo.pinwheel_stats.x(~isnan(PwInfo.pinwheel_stats.x));
%  PWyList = PwInfo.pinwheel_stats.y(~isnan(PwInfo.pinwheel_stats.y));

%  local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PWxList, PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
%  local_pw_dens = local_pw_dens./sum(local_pw_dens(ROI)).*sum(PwInfo.pinwheel_stats.probability,'all')/(sum(ROI,'all')/data_info.pix_per_mm^2);

%  plot_mapAbs(local_pw_dens,'Pinwheel Prob. Density [1/mm^2]',max(local_pw_dens(data_obj.ROI),[],'all'),min(local_pw_dens(data_obj.ROI),[],'all'),data_obj.ROI,gca)
%  hold on; contour(ROI,[1 1],'white','linewidth',2)
%  hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','linewidth',linewidth)
% %         hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
% %         hold on
% %         contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
% %         hold on
% %         contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
% %         hold on
% %         for ii = 1:size(PwInfo.pinwheel_stats.x,1)
% %             text(PwInfo.pinwheel_stats.x(ii,1)+2,PwInfo.pinwheel_stats.y(ii,1),num2str(PwInfo.pinwheel_stats.probability(ii)),'Color',color_contur,'FontSize',5)
% %         end
%  xlim(x_range)
%  ylim(y_range)

%  hold on; set(gca,'view',[rotate rotate])
%  yticks([])
%  xticks([])
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'PwCI.ps'])


%  %% plot Ori CI
%  f1 = figure;
%  plot_mapAbs(NoiseMeasure.CI_ori,'CI ori prefrerence [°]',180,0,ROI,gca)
%  xlim(x_range)
%  ylim(y_range)
%  hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','linewidth',linewidth)
%  hold on; set(gca,'view',[rotate rotate])
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'OriCI.ps'])
 

%  %% plot dir map big
%  f1=figure;
%  plot_map(z_dir,data_obj_dir.ROI,0,1)
%  [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
%  U=cos(angle(z_dir));V=sin(angle(z_dir));
%  U(~ROI) = 0;
%  V(~ROI) = 0;
%  d=10;
%  hold on
%  quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
%  xlim(x_range)
%  ylim(y_range)
%  title('direction map ')
%  hold on; set(gca,'view',[rotate rotate])
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'Dir.ps'])

%  %% plot dir map Abs big
%  f1=figure;
%  plot_Abs(abs(z_dir).*ROI,gca)
%  hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','linewidth',linewidth)
%  hold on
%  quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)

%  hold on
%  contour(real(z_dir),[0 0],'linewidth',linewidth,'Color','white')
%  hold on
%  contour(imag(z_dir),[0 0],'linewidth',linewidth,'Color','white')

%  xlim(x_range)
%  ylim(y_range)
%  title('selectivity direction map')
%  hold on; set(gca,'view',[rotate rotate])
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirAbs.ps'])

%  %% plot Dir CI map
%  f1 = figure;
%  plot_mapAbs(NoiseMeasure.CI_dir,'CI dir prefrerence [°]',360,0,ROI,gca)
%  xlim(x_range)
%  ylim(y_range)
%  hold on; set(gca,'view',[rotate rotate])
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirCI.ps'])

%  %% compare dir selectivity histogramms with norm pdf
%  f1 = figure;
%  histogram(abs(z_dir(ROI)),'Normalization','pdf','NumBins',20,'DisplayName','all sel.')
%  hold on
%  sel_pw = abs(z_dir(round(PwInfo.PWxList),round(PwInfo.PWyList)));
%  histogram(sel_pw,'Normalization','pdf','NumBins',20,'DisplayName','pinwheel sel.')
%  legend()
%  print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirSelHist.ps'])




%         

% 
%         %% make Figure
%         f = figure(1);
%         t = tiledlayout(4,3);
%         title(t,data_info.ID)
            
        
%         %% plot dir map
%         nexttile(t)
%         plot_map(z_dir,data_obj.ROI,0,1)
%         [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
%         U=cos(angle(z_dir));V=sin(angle(z_dir));
%         U(~ROI) = 0;
%         V(~ROI) = 0;
% 
%         d=10;
%         hold on
%         quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
%         xlim([Xmin Xmax])
%         ylim([Ymin Ymax])
%         title('direction map')


%         %% plot dir map data obj
%         z_dir = data_obj_dir.filter_map(data_obj_dir.read_map(1,false,true));
%         nexttile(t)
%         plot_map(z_dir,data_obj_dir.ROI,0,1)
%         [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
%         U=cos(angle(z_dir));V=sin(angle(z_dir));
%         U(~ROI) = 0;
%         V(~ROI) = 0;
% 
%         d=10;
%         hold on
%         quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
%         xlim([Xmin Xmax])
%         ylim([Ymin Ymax])
%         title('direction map ')
% 
%         %% plot ori map
%         nexttile(t)
%         plot_map(z,data_obj.ROI,0,1)
%         [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
%         U=cos(angle(z_dir));V=sin(angle(z_dir));
%         U(~ROI) = 0;
%         V(~ROI) = 0;
% 
%         d=10;
%         hold on
%         quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
%         xlim([Xmin Xmax])
%         ylim([Ymin Ymax])
%         title('orientation map')
% 
%         %% plot ori map data obj
%         z_av = data_obj_dir.filter_map(data_obj_dir.read_map());
%         nexttile(t)
%         plot_map(z_av,data_obj.ROI,0,1)
%         [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
%         U=cos(angle(z_dir));V=sin(angle(z_dir));
%         U(~ROI) = 0;
%         V(~ROI) = 0;
% 
%         d=10;
%         hold on
%         quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
%         xlim([Xmin Xmax])
%         ylim([Ymin Ymax])
%         title('orientation map av dir')
% 
%         [PwInfo_av.NumberPw,PwInfo_av.aniso,PwInfo_av.x_angle,PwInfo_av.PWxList,PwInfo_av.PWyList,PwInfo_av.signList, PwInfo_av.contours] = find_pinwheels(z_av,0,data_obj_dir.ROI);


%         %% plot map average over oppsite directions
%         %DirectionData = DirectionData - mean(DirectionData(:,:,1:(end-1),:),3:4);
%         StimDir = StimDir(:,:,1:end/2,:);
%         z_av = mean((DirectionData(:,:,1:length(StimDir),:)+DirectionData(:,:,length(StimDir)+1:size(DirectionData,3)-1,:))/2.*exp(2i*pi*StimDir/180),3:4);
%         z_av = data_obj.filter_map(z_av);
%         [PwInfo_av.NumberPw,PwInfo_av.aniso,PwInfo_av.x_angle,PwInfo_av.PWxList,PwInfo_av.PWyList,PwInfo_av.signList, PwInfo_av.contours] = find_pinwheels(z_av,0,data_obj.ROI);
% 
%         nexttile(t)
%         plot_map(z_av,data_obj.ROI,0,1)
%         [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
%         U=cos(angle(z_dir));V=sin(angle(z_dir));
%         U(~ROI) = 0;
%         V(~ROI) = 0;
% 
%         d=10;
%         hold on
%         quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
%         xlim([Xmin Xmax])
%         ylim([Ymin Ymax])
%         title('averaged direction map')


%         %% plot CI direction map
%         ax = nexttile(t);
%         plot_mapAbs(NoiseMeasure.CI_dir,'CI dir prefrerence [°]',360,0,ROI,ax)
%         xlim(x_range)
%         ylim(y_range)
% 
% 
% 
% 
% 
%          %% plot CI ori map
%         [NoiseMeasure.CI_ori,~,~] = getCI(data_obj,alpha,'bca',true,false);
%         ax = nexttile(t);
%         plot_mapAbs(NoiseMeasure.CI_ori,'CI ori prefrerence [°]',180,0,ROI,ax)
%         xlim(x_range)
%         ylim(y_range)
% 
%         %% plot CI ori map
%         [NoiseMeasure.CI_ori_av,~,~] = getCI(data_obj_dir,alpha,'bca',true,false);
%         ax = nexttile(t);
%         plot_mapAbs(NoiseMeasure.CI_ori_av,'CI ori prefrerence av [°]',180,0,ROI,ax)
%         xlim(x_range)
%         ylim(y_range)
% 
%         %% plot histogramm CI_dir
%         nexttile(t)
%         histogram(NoiseMeasure.CI_dir(ROI),100)
%         title('CI dir prefrerence [°]')
%         xlim([0 360])
% 
%         %% plot histogramm CI_ori
%         nexttile(t)
%         histogram(NoiseMeasure.CI_ori(ROI),100)
%         title('CI ori prefrerence [°]')
%         xlim([0 180])
% 
%         %% plot histogramm CI_ori_av
%         nexttile(t)
%         histogram(NoiseMeasure.CI_ori_av(ROI),100)
%         title('CI ori prefrerence av [°]')
%         xlim([0 180])
% 
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%% set alpha for SE -> 68.2 CI
%         alpha = 1-0.682;
% 
%                 %% plot CI direction map
%         [x_range,y_range] = getRangeXY_ROI(data_obj_dir.ROI);
%         [NoiseMeasure_SE.CI_dir,~,~] = getCI(data_obj_dir,alpha,'bca',true,true);
%         ax = nexttile(t);
%         plot_mapAbs(NoiseMeasure_SE.CI_dir,'SE dir prefrerence [°]',360,0,ROI,ax)
%         xlim(x_range)
%         ylim(y_range)
% 
% 
%          %% plot CI ori map
%         [NoiseMeasure_SE.CI_ori,~,~] = getCI(data_obj,alpha,'bca',true,false);
%         ax = nexttile(t);
%         plot_mapAbs(NoiseMeasure_SE.CI_ori,'SE ori prefrerence [°]',180,0,ROI,ax)
%         xlim(x_range)
%         ylim(y_range)
% 
%         %% plot CI ori map
%         [NoiseMeasure_SE.CI_ori_av,~,~] = getCI(data_obj_dir,alpha,'bca',true,false);
%         ax = nexttile(t);
%         plot_mapAbs(NoiseMeasure_SE.CI_ori_av,'SE ori prefrerence av [°]',180,0,ROI,ax)
%         xlim(x_range)
%         ylim(y_range)
% 
%         %% plot histogramm CI_dir
%         nexttile(t)
%         histogram(NoiseMeasure_SE.CI_dir(ROI),100)
%         title('SE dir prefrerence [°]')
%         xlim([0 360])
% 
%         %% plot histogramm CI_ori
%         nexttile(t)
%         histogram(NoiseMeasure_SE.CI_ori(ROI),100)
%         title('SE ori prefrerence [°]')
%         xlim([0 180])
% 
%         %% plot histogramm CI_ori_av
%         nexttile(t)
%         histogram(NoiseMeasure_SE.CI_ori_av(ROI),100)
%         title('SE ori prefrerence av [°]')
%         xlim([0 180])
% 
% 




function ResponseData = getResponseData(data,Filter,data_obj)
    ResponseData = zeros(size(data,[1 2 3]));
    for ii = 1:size(ResponseData,3)
        if Filter
            ResponseData(:,:,ii) = data_obj.filter_map(mean(data(:,:,ii,:),4));
        else
            ResponseData(:,:,ii) = mean(data(:,:,ii,:),4);
        end
    end
end

function plot_Abs(map,ax)

    a = map;   
    
 
    imagesc(ax,a);

    %% details Plot
    %colormap jet;
    %colormap turbo;
    m=100;
    cm_viridis=viridis(m);
    colormap(ax, cm_viridis);%gray,'jet'
    %colormap hot;
    
    set(gca,'xtick',[])
    set(gca,'ytick',[])
%     title(ax,Title)%, 'interpreter', 'latex'
    %pbaspect([1 1 1])
    
    hold on
    %colorbar('Ticks',[0,0.5,3])
    cbh = colorbar(ax) ; %Create Colorbar
%     axis equal
    axis image
%     cbh.Ticks = linspace(0, 1, NTicks) ; %Create 8 ticks from zero to 1
%     cbh.TickLabels = num2cell(linspace(minMap,maxMap,NTicks)) ;
    %cbh.TickLabels = num2cell(linspace(-1,1,NTicks)) ;
end
    
    


function LocalPwDensityFixedFilter = getLocalPwDensityFixedFilter(data_obj,PwInfo,local_spacing_mm,sigma)
    disp('calc LocalPwDensity with FixedFilter')
    average_spacing_mm = mean(local_spacing_mm(data_obj.ROI));
    local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PwInfo.PWxList, PwInfo.PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
    local_pw_dens = local_pw_dens./sum(local_pw_dens(data_obj.ROI)).*PwInfo.NumberPw;
    LocalPwDensityFixedFilter = local_pw_dens.*(local_spacing_mm*data_obj.info.pix_per_mm).^2;
end

function overlay_white_shade(A)
% OVERLAY_WHITE_SHADE Overlays a white transparent shade on the current plot
% The opacity of the overlay is inversely proportional to the values in A.
%
% Input:
%   A - 2D array with numeric values (does not need to be normalized)
%
% Usage:
%   imagesc(peaks(100));
%   overlay_white_shade(A);

    % Normalize A to range [0, 1]
    A = double(A);  % ensure it's double for math operations
    A = A - min(A(:));
    if max(A(:)) > 0
        A = A ./ max(A(:));
    end
    
    % Compute alpha channel (inverse of A)
    alpha_mask = 1 - A;

    % Prepare white RGB image
    white_overlay = ones([size(A), 3]);

    % Hold current plot and overlay
    hold on;
    h = imagesc(white_overlay);
    set(h, 'AlphaData', alpha_mask);
    uistack(h, 'top');  % Ensure it's on top of existing content
    hold off;
end