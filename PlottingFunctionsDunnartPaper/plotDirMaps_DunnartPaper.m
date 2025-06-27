function plotDirMaps_DunnartPaper(experiments,AnimalDataFolder,DataFolderMain,...
    GIF_SN_TH,data_info_file_list,smallest_w_mm,largest_w_mm,w_step_mm,Bootstrapsamples, plotComparison)
    % experiments.dunnart = [4]
    % plotDirMaps_DunnartPaper(experiments,'~/CIDBN/','/home/michael/','GIF',4,0.2,2,0.1,100)
    
%     experiments.cat_jung = 1:5;
%     experiments.microcebus = 1:4;
%     experiments.dunnart = [4 1 2 5 6];
%     plotDirMaps_DunnartPaper(experiments,'~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/OrientationPrefernceMapProcessing/AllData&Results/DirectionDataPlots/',4,'experiment_info.csv',0.1,2,0.05,100)

    close all
    %% parameter
    %experiment_num_list = checkFormatNum(experiment_num_list);
    if nargin <4
        GIF_SN_TH = 4;
    end
    if nargin < 5 %|| isempty(data_info_file_list)
        function_dir = fileparts(mfilename('fullpath')); % Get function directory
        data_info_file_list = fullfile(function_dir, 'experiment_info.csv');
    end

    if nargin < 10
        plotComparison = true;
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
    result_data_coll = struct();

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
            result_data.(experiment){experiment_num} = plotOneAnimal(experiment,experiment_num,AnimalDataFolder,DataFolder,data_info_file_list,GIF_SN_TH,smallest_w_mm,largest_w_mm,w_step_mm,Bootstrapsamples,alpha);
        end

        if plotComparison
            %% comparison plots one species
            result_data_coll.(experiment) = struct();
            %% SE distributions CPDF
            f1 = figure;
            f1.Position = [1 1 200 200];
            result_data_coll.(experiment).CI_dir = [];
            for jj = 1:length(experiment_num_list)
                experiment_num = experiment_num_list(jj);
                data = result_data.(experiment){experiment_num};
                CI_dir = data.data_dir.CI;
                ROI = data.data_dir.ROI;
                CI_dir = CI_dir(ROI);
                plot(sort(CI_dir),linspace(0,1,length(CI_dir)),'DisplayName',[experiment num2str(experiment_num)])
                hold on
                result_data_coll.(experiment).CI_dir = [result_data_coll.(experiment).CI_dir; CI_dir];
            end
            ylabel('cum. prob.')
            xlabel(['SE direction sel. [°]'])
            xlim([0 360])
            print(f1,'-dpsc', [DataFolderMain experiment 'SE_CDF.ps'])
    
            %% violine plots modularity distributions Dir
            peak_test_list = cell(1,length(experiment_num_list));
            peak_test_rand_list = cell(1,length(experiment_num_list));
            result_data_coll.(experiment).peaks_test = [];
            result_data_coll.(experiment).peaks_test_rand = [];
            for jj = 1:length(experiment_num_list)
                experiment_num = experiment_num_list(jj);
                data = result_data.(experiment){experiment_num};
                peak_test_list{jj} = permute(data.peaks_test,[2 1])/max(data.peaks_test);
                peak_test_rand_list{jj} = permute(data.peaks_test_rand,[2 1])/max(data.peaks_test);
    
                result_data_coll.(experiment).peaks_test = [result_data_coll.(experiment).peaks_test data.peaks_test/max(data.peaks_test)];
                result_data_coll.(experiment).peaks_test_rand = [result_data_coll.(experiment).peaks_test_rand data.peaks_test_rand/max(data.peaks_test)];
            end
    
            
            f1 = figure;
            f1.Position = [1 1 350/2*length(peak_test_list) 350];
            % plot peak test in blue
            violin(peak_test_list,'facecolor',[0 0 1]);
            % plot rand test in red
            hold on
            violin(peak_test_rand_list,'facecolor',[1 0 0]);
            ylim([0 1.5])
            % set x- ticks to 1,2,3...
            xticks(1:length(experiment_num_list))
            xticklabels(experiment_num_list)
            xlabel(experiment)
            ylabel('Modularity Peak Hight')
            print(f1,'-dpsc', [DataFolderMain experiment 'Modularity_Violine.ps'])
        
            if (experiment == "dunnart")
                %% column spacing dir vs ori
                % get column spacing CI
                upper_dir = [];
                lower_dir = [];
                mean_dir = [];
                upper_ori = [];
                lower_ori = [];
                mean_ori = [];
                for jj = 1:length(experiment_num_list)
                    experiment_num = experiment_num_list(jj);
                    data = result_data.(experiment){experiment_num};
                    upper_dir = [upper_dir; data.data_dir.CI_average_spacing_mm(2)];
                    lower_dir = [lower_dir; data.data_dir.CI_average_spacing_mm(1)];
                    mean_dir = [mean_dir; data.data_dir.average_spacing_mm];
                    upper_ori = [upper_ori; data.data_ori.CI_average_spacing_mm(2)];
                    lower_ori = [lower_ori; data.data_ori.CI_average_spacing_mm(1)];
                    mean_ori = [mean_ori; data.data_ori.average_spacing_mm];
                end
                f1 = figure;
                f1.Position = [1 1 200 200];
                % plot dir y-axis and ori x-axis with error bars
                errorbar(mean_ori,mean_dir,mean_dir-lower_dir,upper_dir-mean_dir,mean_ori-lower_ori,upper_ori-mean_ori,'o')
                xlabel('Ori Column Spacing [mm]')
                ylabel('Dir Column Spacing [mm]')
                print(f1,'-dpsc', [DataFolderMain experiment 'ColumnSpacingDirOri.ps'])
            end
            %% plot angle ori vs dir in degree
            dir_angle_list = [];
            ori_angle_list = [];
            for jj = 1:length(experiment_num_list)
                experiment_num = experiment_num_list(jj);
                data = result_data.(experiment){experiment_num};
                CI_dir = data.data_dir.CI;
    %             dir_angle_list = [dir_angle_list; angle(data.data_dir.z_dir(data.data_dir.ROI))];
    %             ori_angle_list = [ori_angle_list; angle(data.data_ori.z(data.data_ori.ROI))];
                indx = CI_dir<90;
                %indx = abs(data.data_dir.z_dir)<mean(abs(data.data_dir.z_dir(data.data_ori.ROI)))*1.1;
                dir_angle_list = [dir_angle_list; angle(data.data_dir.z_dir(indx))];
                ori_angle_list = [ori_angle_list; angle(data.data_ori.z(indx))];
            end
    
            % convert to degree
            dir_angle_list = dir_angle_list/2/pi*360;
            ori_angle_list = ori_angle_list/2/pi*180;
    
            % scatter plot
            f1 = figure;
            f1.Position = [1 1 200 200];
            plot(ori_angle_list,dir_angle_list,'.','MarkerSize',0.05)
            
    %         xedges = linspace(-90, 90, 50);
    %         yedges = linspace(-180, 180, 50);
    %         
    %         % Compute 2D histogram
    %         N = histcounts2(ori_angle_list, dir_angle_list, xedges, yedges);
    %         plot_mapAbs(N,'angle distributions',max(N(:)),min(N(:)),ones(size(N)),gca)
            
            xlim([-90 90])
            ylim([-180 180])
            xlabel('Ori Angle [°]')
            ylabel('Dir Angle [°]')
            print(f1,'-dpsc', [DataFolderMain experiment 'AngleDirOri.ps'])
    
            %% plot histogramm difference angles ori and dir
            f1 = figure;
            f1.Position = [1 1 200 200];
            histogram(mod(dir_angle_list-ori_angle_list,180),'Normalization','pdf','NumBins',20,'DisplayName','all sel.')
            xlabel('angle difference [°]')
            ylabel('PDF')
            title('angle difference dir - ori')
            print(f1,'-dpsc', [DataFolderMain experiment 'AngleDiffDirOri.ps'])
    
        end

        %% collect data for selectivity vs pw selectivit
        result_data_coll.(experiment).dir_abs = [];
        result_data_coll.(experiment).dir_abs_Pw = [];
        for jj = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(jj);
            data = result_data.(experiment){experiment_num};
            result_data_coll.(experiment).dir_abs = abs(data.data_dir.z_dir(data.data_dir.ROI));
            result_data_coll.(experiment).dir_abs_Pw = abs(data.data_dir.z_dir(sub2ind(size(data.data_dir.z_dir),round(data.data_ori.PwInfo.PWyList),round(data.data_ori.PwInfo.PWxList))));
        end

        %% collect dir angle data
        result_data_coll.(experiment).dir_angle = [];
        for jj = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(jj);
            data = result_data.(experiment){experiment_num};
            result_data_coll.(experiment).dir_angle = [result_data_coll.(experiment).dir_angle; angle(data.data_dir.z_dir(data.data_dir.ROI))];
        end

        %% collect dir gradient data
        result_data_coll.(experiment).inters_angles = [];
        result_data_coll.(experiment).angle_weights = [];

        for jj = 1:length(experiment_num_list)
            experiment_num = experiment_num_list(jj);
            data = result_data.(experiment){experiment_num};
            result_data_coll.(experiment).inters_angles = [result_data_coll.(experiment).inters_angles; data.inters_angles];
            result_data_coll.(experiment).angle_weights = [result_data_coll.(experiment).angle_weights; data.angle_weights];
        end
        
    
    end

    if plotComparison
        %% violine plot modularity all animals
        peak_test_list = cell(1,length(experiment_list));
        peak_test_rand_list = cell(1,length(experiment_list));
        xlabel_list = cell(1,length(experiment_list));
        for jj = 1:length(experiment_list)
            experiment = experiment_list{jj};
            peak_test_list{jj} = permute(result_data_coll.(experiment).peaks_test,[2 1]);
            peak_test_rand_list{jj} = permute(result_data_coll.(experiment).peaks_test_rand,[2 1]);
            xlabel_list{jj} = experiment;
        end
    
        f1 = figure;
        % plot peak test in blue
        violin(peak_test_list,'facecolor',[0 0 1],'xlabel',xlabel_list);
        % plot rand test in red
        hold on
        violin(peak_test_rand_list,'facecolor',[1 0 0],'xlabel',xlabel_list);
        ylim([0 1.5])
        print(f1,'-dpsc', [DataFolderMain 'Modularity_Violine_All.ps'])
    
    
    
    
        %% SE distributions CPDF
        f1 = figure;
        f1.Position = [1 1 200 200];
        for ii = 1:length(experiment_list)
            experiment = experiment_list{ii};
            CI_dir = result_data_coll.(experiment).CI_dir;
            plot(sort(CI_dir),linspace(0,1,length(CI_dir)),'DisplayName',experiment)
            hold on
        end     
        legend()
        xlabel('SE')
        ylabel('CPDF')
        print(f1,'-dpsc', [DataFolderMain 'SE_CDF_All.ps'])
    
        %% compare dir selectivity histogramms pw with norm pdf
        f1 = figure;
        f1.Position = [1 1 200 200];
        %set(f1, 'InvertHardcopy', 'off'); 
        % histogram(abs(z_dir(ROI)),'Normalization','pdf','NumBins',20,'DisplayName','all sel.')
        % hold on
        % sel_pw = abs(z_dir(round(PwInfo.PWxList),round(PwInfo.PWyList)));
        % histogram(sel_pw,'Normalization','pdf','NumBins',20,'DisplayName','pinwheel sel.')
        % legend()
        % print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirSelHist.ps'])
        for ii = 1:length(experiment_list)
            experiment = experiment_list{ii};
            data = result_data_coll.(experiment);
            % histogram(data.dir_abs/mean(data.dir_abs),'Normalization','pdf','NumBins',20,'DisplayName',[experiment ' all sel.'],'DisplayStyle', 'stairs')
            % hold on
            sel_pw = data.dir_abs_Pw;
            histogram(sel_pw/mean(data.dir_abs),'Normalization','pdf','NumBins',10,'DisplayName',[experiment ' pinwheel sel.'],'DisplayStyle', 'stairs')%,'FaceAlpha',0.3
            hold on
        end
        legend()
        xlabel('rel. Dir Selectivity')
        ylabel('PDF')
        print(f1,'-dpsc', [DataFolderMain 'DirSelHistPw_All.ps']), %'-image'

        %% compare dir selectivity histogramms with norm pdf
        f1 = figure;
        f1.Position = [1 1 200 200];
        %set(f1, 'InvertHardcopy', 'off'); 
        % histogram(abs(z_dir(ROI)),'Normalization','pdf','NumBins',20,'DisplayName','all sel.')
        % hold on
        % sel_pw = abs(z_dir(round(PwInfo.PWxList),round(PwInfo.PWyList)));
        % histogram(sel_pw,'Normalization','pdf','NumBins',20,'DisplayName','pinwheel sel.')
        % legend()
        % print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirSelHist.ps'])
        for ii = 1:length(experiment_list)
            experiment = experiment_list{ii};
            data = result_data_coll.(experiment);
            histogram(data.dir_abs/mean(data.dir_abs),'Normalization','pdf','NumBins',20,'DisplayName',[experiment ' all sel.'],'DisplayStyle', 'stairs')
            hold on
            % sel_pw = data.dir_abs_Pw;
            % histogram(sel_pw/mean(data.dir_abs),'Normalization','pdf','NumBins',20,'DisplayName',[experiment ' pinwheel sel.'],'DisplayStyle', 'stairs')%,'FaceAlpha',0.3
        end
        legend()
        xlabel('rel. Dir Selectivity')
        ylabel('PDF')
        xlim([0 5])
        ylim([0 1])
        print(f1,'-dpsc', [DataFolderMain 'DirSelHist_All.ps']), %'-image'
    
        %% compare angle data dir
        f1 = figure;
        f1.Position = [1 1 300 300];
        for ii = 1:length(experiment_list)
            experiment = experiment_list{ii};
            data = result_data_coll.(experiment);
            histogram(rad2deg(data.dir_angle),'Normalization','pdf','NumBins',20,'DisplayName',[experiment],'DisplayStyle', 'stairs')
            hold on
    
        end
        xlabel('Dir Angle [°]')
        ylabel('PDF')
        xlim([0 5])
        ylim([0 1])
        legend()
        print(f1,'-dpsc', [DataFolderMain 'DirAngleHist_All.ps'])


        %% compare inters angles with weighted histograms
        f1 = figure;
        f1.Position = [1 1 300 300];
        for ii = 1:length(experiment_list)
            experiment = experiment_list{ii};
            data = result_data_coll.(experiment);
            inters_angles = data.inters_angles/2/pi*360;
            wheigths = data.angle_weights/sum(data.angle_weights);
            weightedhistogram(inters_angles,wheigths,'Label',experiment, 'Normalize', 'pdf')
            hold on
        end
        legend()
        xlabel('intersection angle [°]')
        ylabel('PDF')
        print(f1,'-dpsc', [DataFolderMain 'IntersectionAngles.ps'])

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

    %% define ResultDataFile
    ResultDataFile = [DataFolder animal num2str(experiment_num) 'ResultData.mat'];
    
    if isfile(ResultDataFile)
        disp(['ResultDataFile already exists: ' ResultDataFile])
        load(ResultDataFile,'result_data','data_info');
        disp('loading data from file')
    else
   
        %% get DirData
        [result_data.data_dir,data_obj_dir] = getDirResultData(data_info,data_path,data_obj,GIF_SN_TH,Bootstrapsamples,alpha);

        %% get OriData
        [result_data.data_ori] = getOriResultData(data_obj,alpha);


        %% get CI column spacing Dir
        if animal == "dunnart"
            getCI = true;
            %[average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap,alpha)
            [result_data.data_dir.average_spacing_mm,result_data.data_dir.local_spacing_mm,~,~,result_data.data_dir.CI_average_spacing_mm,result_data.data_dir.CI_local_spacing_mm] = getColumnsSpacing(data_obj_dir,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,true,alpha,true);

            %% get CI column spacing Ori
            %[average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap,alpha)   
            [result_data.data_ori.average_spacing_mm,result_data.data_ori.local_spacing_mm,~,~,result_data.data_ori.CI_average_spacing_mm,result_data.data_ori.CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,true,alpha,false);

        
        else
            getCI = false;
            %[average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap,alpha)
            [result_data.data_dir.average_spacing_mm,result_data.data_dir.local_spacing_mm,~,~,~,~] = getColumnsSpacing(data_obj_dir,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,true,alpha,true);

            %% get CI column spacing Ori
            %[average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap,alpha)   
            [result_data.data_ori.average_spacing_mm,result_data.data_ori.local_spacing_mm,~,~,~,~] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,true,alpha,false);

        
        end


        %% test modularity
        profile_range_mm = smallest_w_mm:w_step_mm:largest_w_mm;
        [result_data.peaks_test,result_data.peaks_test_rand,result_data.peak_position_mm,result_data.power_profiles,result_data.power_profiles_rand,~,~] = testModularityOPM(data_obj_dir,DataFolder,result_data.data_dir.average_spacing_mm,profile_range_mm,Bootstrapsamples,false,true);


        %% get iso-line intersetcion angles and gradients
        [result_data.grad_dir,result_data.grad_ori,result_data.inters_angles,result_data.angle_weights] = getIsoLineAngles(result_data.data_dir.z_dir,result_data.data_ori.z,data_obj.ROI);

        %% save data
        save(ResultDataFile,'result_data','data_obj_dir','data_info','-v7.3');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% make map borders ROI
    rel_x_length = 2;
    rel_y_length = 3;
    [x_range,y_range] = makeRectangleROI(data_obj,result_data.data_ori.average_spacing_mm,rel_x_length,rel_y_length);

    


    %% set all variables
    z = result_data.data_ori.z;
    z_dir = result_data.data_dir.z_dir;
    ROI = data_obj.ROI;
    CI_dir = result_data.data_dir.CI;
    PwInfo = result_data.data_ori.PwInfo;


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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% plot dir polar map + arrows + scale bar 0.5 mm
    f1 = figure();
    plot_map_DIR(z_dir,ROI,3*sqrt(mean(abs(z_dir(ROI)).^2)),1)
    xlim(x_range)
    ylim(y_range)
    title('direction pref. polar map ')

    [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
    U=cos(angle(z_dir));V=sin(angle(z_dir));
    U(~ROI) = 0;
    V(~ROI) = 0;
    d=10;
    hold on
    quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)

    hold on
    Xmin = x_range(1);
    Ymin = y_range(1);
    mm = 0.5;
    width_scale_pix = 5;
    spacing_pix = mm * data_info.pix_per_mm;
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-white','LineWidth', 1.5)

    if ~isempty(ElectrodePosition)
        hold on;
        plot(ElectrodePosition(1),ElectrodePosition(2),'xwhite','MarkerSize',15)
    end

    hold on; set(gca,'view',[rotate rotate])
    print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'PolarDir.ps'])

    %% plot dir polar map + arrows + scale bar 0.5 mm
    f1 = figure();
    plot_map_DIR(z_dir,ROI,3*sqrt(mean(abs(z_dir(ROI)).^2)),1)
    xlim(x_range)
    ylim(y_range)
    title('direction pref. polar map ')

    [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
    U=cos(angle(z_dir));V=sin(angle(z_dir));
    U(~ROI) = 0;
    V(~ROI) = 0;
    d=10;
    hold on
    quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)

    hold on
    Xmin = x_range(1);
    Ymin = y_range(1);
    mm = 0.5;
    width_scale_pix = 5;
    spacing_pix = mm * data_info.pix_per_mm;
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-white','LineWidth', 1.5)

    if ~isempty(ElectrodePosition)
        hold on;
        plot(ElectrodePosition(1),ElectrodePosition(2),'xwhite','MarkerSize',15)
    end

    hold on; set(gca,'view',[rotate rotate])
    print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'PolarDir.ps'])

    %% plot angle map Dir + Dir CI
    f1=figure;
    plot_map_DIR(z_dir,ROI,0,1)

    hold on; contour(CI_dir,[80 80],'white','linewidth',2,'LineStyle','-')

    hold on; contour(CI_dir,[40 40],'white','linewidth',2,'LineStyle',':')


    hold on
    Xmin = x_range(1);
    Ymin = y_range(1);
    mm = result_data.data_dir.average_spacing_mm;
    width_scale_pix = 5;
    spacing_pix = mm * data_info.pix_per_mm;
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-white','LineWidth', 1.5)

    xlim(x_range)
    ylim(y_range)
    title(['dirmap + SE + CS= ' num2str(round(result_data.data_dir.average_spacing_mm,3)) ' mm'])
    hold on; set(gca,'view',[rotate rotate])
    print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirCI.ps'])

    %% plot dir polar map + CI
    f1 = figure();
    plot_map_DIR(z_dir,ROI,3*sqrt(mean(abs(z_dir(ROI)).^2)),1)
    xlim(x_range)
    ylim(y_range)
    title('direction pref. polar map SE ')

    hold on; contour(CI_dir,[60 60],'white','linewidth',2,'LineStyle','-')

    hold on; set(gca,'view',[rotate rotate])
    print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'PolarDirCI.ps'])


    %% plot zero comntors ori and dir map on black background
    f1 = figure();

    % background in black
    imagesc(zeros(size(z_dir)),'AlphaData',0.5);
    
    contour(real(z_dir),[0 0],'black','linewidth',2,'LineStyle','--')
    hold on; contour(imag(z_dir),[0 0],'black','linewidth',2,'LineStyle','--')
    
    hold on; contour(real(z),[0 0],'black','linewidth',2,'LineStyle','-')
    hold on; contour(imag(z),[0 0],'black','linewidth',2,'LineStyle','-')

    hold on; contour(ROI,[1 1],'red','linewidth',2,'LineStyle','-')

    xlim(x_range)
    ylim(y_range)
    hold on; set(gca,'view',[rotate rotate])

    print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'ZeroContours.ps'])

    

    

    %% plot dir CI

    f1=figure;

    cb = plot_mapAbs(CI_dir,'selectivity direction map',max(CI_dir(:)),min(CI_dir(:)),ROI,gca);
    
    xlim(x_range)
    ylim(y_range)
    title('SE direction prefernce map ')
    hold on; set(gca,'view',[rotate rotate])

    ax = gca;
    ax_pos = ax.Position;        % [x y width height]
    cb.Position(2) = ax_pos(2);  % align vertical position
    cb.Position(4) = ax_pos(4); 

    print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirCIPref.ps'])

    %% angle map Ori + arrows Dir + Pw positions
    f1 = figure();
    plot_map(z,ROI,0,1)
    xlim(x_range)
    ylim(y_range)
    title('ori pref. + arrows dir map ')

    [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
    U=cos(angle(z_dir));V=sin(angle(z_dir));
    U(~ROI) = 0;
    V(~ROI) = 0;
    d=10;
    hold on
    quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 1)

    hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','MarkerSize',20)

    hold on; set(gca,'view',[rotate rotate])
    print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'OriDir.ps'])

    %% angle map Ori +fracture Lines
    f1 = figure();           % white background
    set(f1, 'InvertHardcopy', 'off'); 


    plot_map(z,ROI,0,1)
    xlim(x_range)
    ylim(y_range)
    title('ori pref. + arrows dir map ')

    [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
    U=cos(angle(z_dir));V=sin(angle(z_dir));
    U(~ROI) = 0;
    V(~ROI) = 0;
    d=5;
    hold on
    quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 1)

    hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','MarkerSize',15)

    overlay_white_shade(abs(z_dir),.9,35,ROI)

    hold on; set(gca,'view',[rotate rotate])
%     print(f1,'-dpsc','-fillpage', '-image', [DataFolder animal num2str(experiment_num) 'OriDirFract.ps'])
    print(f1, '-dpsc', '-image', '-fillpage', '-r1200', ...
    [DataFolder animal num2str(experiment_num) 'OriDirFract.ps']);
    


%     %% plot Abs Dir map + Pws
%     f1=figure;
%     f1.Position = [1 1 400 400];
%     cb = plot_mapAbs(abs(z_dir),'selectivity direction map',max(abs(z_dir(:))),0,ROI,gca);
% 
%     xlim(x_range)
%     ylim(y_range)
% 
%     [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
%     U=cos(angle(z_dir));V=sin(angle(z_dir));
%     U(~ROI) = 0;
%     V(~ROI) = 0;
%     d=10;
%     hold on
%     quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
% 
% 
%     hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','MarkerSize',10)
% 
%     title('selectivity direction map')
%     
% 
%     ax = gca;
%     ax_pos = ax.Position;        % [x y width height]
%     cb.Position(2) = ax_pos(2);  % align vertical position
%     cb.Position(4) = ax_pos(4); 
% 
%     hold on; set(gca,'view',[rotate rotate])
% 
%     print(f1,'-dpsc', [DataFolder animal num2str(experiment_num) 'DirAbs.ps'])
    
    %% plot Abs Dir map + Pws
    f1 = figure();
    
    imagesc(gca,abs(z_dir).*ROI,[0 max(abs(z_dir(:)))]);
    m=100;
    cm_viridis=viridis(m);
    colormap(ax, cm_viridis)
    
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

    [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
    U=cos(angle(z_dir));V=sin(angle(z_dir));
    U(~ROI) = 0;
    V(~ROI) = 0;
    d=10;
    hold on
    quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)


    hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','MarkerSize',10)

    

    ax = gca;
    ax_pos = ax.Position;        % [x y width height]
    cb.Position(2) = ax_pos(2);  % align vertical position
    cb.Position(4) = ax_pos(4); 

    hold on; set(gca,'view',[rotate rotate])
    
    title('selectivity direction map')

    print(f1,'-dpsc', [DataFolder animal num2str(experiment_num) 'DirAbs.ps'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     %% get CI column spacing Dir
%     %[average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap,alpha)
%     [result_data.data_dir.average_spacing_mm,result_data.data_dir.local_spacing_mm,~,~,result_data.data_dir.CI_average_spacing_mm,result_data.data_dir.CI_local_spacing_mm] = getColumnsSpacing(data_obj_dir,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,true,true,alpha,true);
% 
%     %% get CI column spacing Ori
%     %[average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap,alpha)   
%     [result_data.data_ori.average_spacing_mm,result_data.data_ori.local_spacing_mm,~,~,result_data.data_ori.CI_average_spacing_mm,result_data.data_ori.CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,true,true,alpha,false);



end

function [grad_dir,grad_ori,inters_angles,weights] = getIsoLineAngles(z_dir,z_ori,ROI)

    %% get gradients
    dx = 1;
    dy = 1;
    [dir_dtheta_dx, dir_dtheta_dy] = phase_gradient(z_dir, dx, dy);
    [ori_dtheta_dx, ori_dtheta_dy] = phase_gradient(z_ori, dx, dy);

    %% extract values inside ROI
    dir_dtheta_dx = dir_dtheta_dx(ROI);
    dir_dtheta_dy = dir_dtheta_dy(ROI);
    ori_dtheta_dx = ori_dtheta_dx(ROI);
    ori_dtheta_dy = ori_dtheta_dy(ROI);

    %% calc angles
    scalar_product = dir_dtheta_dx .* ori_dtheta_dx + dir_dtheta_dy .* ori_dtheta_dy;
    dir_norm = sqrt(dir_dtheta_dx.^2 + dir_dtheta_dy.^2);
    ori_norm = sqrt(ori_dtheta_dx.^2 + ori_dtheta_dy.^2);
    cos_angle = scalar_product ./ (dir_norm .* ori_norm);
    inters_angles = acos(cos_angle); % angles in radians

    %% set gradients
    grad_dir.dx = dir_dtheta_dx;
    grad_dir.dy = dir_dtheta_dy;
    grad_ori.dx = ori_dtheta_dx;
    grad_ori.dy = ori_dtheta_dy;

    %% set weights
    weights = ori_norm .* dir_norm;

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


function [result_data] = getOriResultData(data_obj,alpha)%smallest_w_mm,largest_w_mm,w_step_mm

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

function [result_data,data_obj_dir]=getDirResultData(data_info,data_path,data_obj,GIF_SN_TH,Bootstrapsamples,alpha)
    if data_info.animal == "Mouse Lemur"
            data_obj_dir = data_obj;
    else
        %% get dir map
        [~,~,DirectionData] = getDirectionData(data_info,data_path,0);

        DirectionData  = real(DirectionData);
    
        %% make data obj dir data
        data_info_dir = data_info;
        stimDir = [data_info.stim_order(find(~isnan(data_info.stim_order))) data_info.stim_order(find(~isnan(data_info.stim_order)))+180];
        stimDir = [stimDir NaN];
        data_info_dir.stim_order=stimDir;
        data_obj_dir = data_handle_corrected(data_info_dir,DirectionData,data_obj.ROI);
    end
    
    %% prepare GIF
    if GIF_SN_TH>0
        data_obj_dir.activateGIF(true,GIF_SN_TH)
    end

    %% prepare bootstrap samples
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

function overlay_white_shade(A,step,k,ROI)
% OVERLAY_WHITE_SHADE Overlays a white transparent shade on the current plot
% The opacity of the overlay is inversely proportional to the values in A.
%
% Input:
%   A - 2D array with numeric values (does not need to be normalized)
%
% Usage:
%   imagesc(peaks(100));
%   overlay_white_shade(A);

%   overlay_white_shade(A, 0.5); % step is optional, default is 0.5

    % Normalize A to range [0, 1]
    A = double(A);  % ensure it's double for math operations
    A = A - min(A(:));
    if max(A(:)) > 0
        A = A ./ max(A(:));
    end
    
    % Compute alpha channel (inverse of A)
    alpha_mask = 1 - A;

    %% apply step function to alpha mask
    
%     alpha_mask(alpha_mask < step) = 0;
%     alpha_mask(alpha_mask >= step) = 1;


    alpha_mask =  smooth_step(alpha_mask, k, step);

    alpha_mask(~ROI)=0;

    % Prepare white RGB image
    white_overlay = ones([size(A), 3]);

    % Hold current plot and overlay
    hold on;
    h = imagesc(white_overlay);
    set(h, 'AlphaData', alpha_mask);
    uistack(h, 'top');  % Ensure it's on top of existing content
    hold off;
end

function y = smooth_step(x, k, step)
    % smooth_step returns a sigmoid-shaped step function
    % Inputs:
    %   x     - input array
    %   k     - steepness (higher = steeper)
    %   step  - optional: location of the step (default = 0)

    if nargin < 3
        step = 0;
    end

    y = 1 ./ (1 + exp(-k * (x - step)));
end

function varargout=plot_map_DIR(data,ROI,ref_sel,black_roi,smoothing,color_type)
% [h,fig2plot]=plot_map(data,ROI,ref_sel,black_roi,smoothing,color_type)
% Function to make a color and intensity colored orientation preference
% display
% input = polar data, Region of interest
% mapPolar(data2D,ROI)
% created by chepe@nld.ds.mpg.de

% Other functions used:
% - makeColormap
% - smooth_pf, negative_smoothing, imresize


%% ---------------------  read parameters

data(isnan(data)) = 0;

% if empty ROI, use all image
if ~exist('ROI','var') || isempty(ROI)
    ROI=true(size(data));
end

% determine polar scaling for selectivity
if ~exist('ref_sel','var')
    ref_sel=3*sqrt(mean(abs(data(ROI)).^2));
elseif isempty(ref_sel)
    ref_sel=3*sqrt(mean(abs(data(ROI)).^2));
end

if ~exist('black_roi','var')
    black_roi = 1;
end

% increase size of display for exporting figures
if exist('smoothing','var') && ~isempty(smoothing)
    if smoothing>0
        for times=1:smoothing
            data=smooth_pf(data);
        end
    elseif smoothing<0
        for times=1:abs(smoothing)
            data = negative_smoothing(data);
        end
    end
    ROI = imresize(ROI,2^abs(smoothing));
end

%-  make color maps
if ~exist('color_type','var')
    color_type = 'std';
end
switch color_type
    case 'std'
        cm = makeColormap2('orientation',16);        
    case 'interp'
        cm = makeColormap2('orientation',24);
    case 'circ'
        cm = makeColormap2('circular',50);        
end
mm = makeColormap2('selectivity',64);
%% --------------------------- Asign color values

% -- ORIENTATION

% normalize from 0-1
ori=(angle(data))/(2*pi);ori(ori<0)=1+ori(ori<0);
ori=ori(:);

% find to which color interval it belongs
intervals=linspace(0,1,2*size(cm,1)+1);
[~,ori] = histc(ori,intervals(2:2:end));
ori = ori+1;
ori(ori==size(cm,1)+1)=1;

% assign
anglePlot = (cm(ori(:),:));

% -- SELECTIVITY

if ref_sel==0
    
    % if polar plot is not required
    magnitudePlot = 0.9*ones(size(anglePlot));
    
else
    
    % threshold
    sel=abs(data)/ref_sel;
    sel(sel>1)=1;
    
    % find intensity match
    sel=ceil(size(mm,1)*sel(:));
    sel(sel==0)=1;
    
    % assign
    magnitudePlot = (mm(sel(:),:));
end

%% --------------------------- Combine

if ismatrix(data)
    magnitudePlot = reshape(magnitudePlot, [size(data,1),size(data,2),3]);
    anglePlot = reshape(anglePlot, [size(data,1),size(data,2),3]);
    fig2plot=anglePlot.*magnitudePlot;
    %fig2plot=repmat(double(ROI),[1 1 3]).*anglePlot.*magnitudePlot;
    if black_roi == 1
        % black ROI
        fig2plot(repmat(~ROI,[1 1 3])) = 0;
    elseif black_roi == -1
        % white ROI
        fig2plot(repmat(~ROI,[1 1 3])) = 1;
    end
else
    % in case of volume image, return vector
    fig2plot = magnitudePlot.*anglePlot;
end

% if only the colormap is required, return without plotting
if nargout==2 || ~ismatrix(data)
    varargout{1} = [];
    varargout{2} = fig2plot;   
    return
end

%% --------------------------- Make figure

h = image(fig2plot);
set(gca,'ydir','reverse')
axis image

set(gca,'xtick',[])
set(gca,'ytick',[])

set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

if nargout==1
    varargout{1} = h;    
end

end
