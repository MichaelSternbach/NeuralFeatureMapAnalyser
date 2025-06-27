function plotOriDirMaps(experiment_num_list,animal,AnimalDataFolder,DataFolderMain,...
    GIF_SN_TH,data_info_file_list)
    
    close all
    %% parameter
    experiment_num_list = checkFormatNum(experiment_num_list);
    if nargin <5
        GIF_SN_TH = 4;
    end
    if nargin < 6 %|| isempty(data_info_file_list)
        function_dir = fileparts(mfilename('fullpath')); % Get function directory
        data_info_file_list = fullfile(function_dir, 'experiment_info.csv');
    end

    rotate = -90;
    linewidth = 1;
    Bootstrapsamples = 100;
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
  

    for experiment_num = experiment_num_list
        
        mkdir([DataFolderMain animal])
        DataFolder = [DataFolderMain animal '/' animal num2str(experiment_num) '/'];
        mkdir(DataFolder)
        
        %% get data
        [data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,experiment_num,AnimalDataFolder,data_info_file_list);
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
        z = data_obj.filter_map(data_obj.read_map());



        %% get dir map
        [DirectionMap,StimDir,DirectionData] = getDirectionData(data_info,data_path,0);
%         z_dir = DirectionMap;
%         z_dir = data_obj.filter_map(DirectionMap);
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
        z_dir = data_obj_dir.filter_map(data_obj_dir.read_map(1,false,true));
    
        %% calc CI dir and Ori map
        [NoiseMeasure.CI_dir,~,~] = getCI(data_obj_dir,alpha,'bca',true,true);
        [NoiseMeasure.CI_ori,~,~] = getCI(data_obj,alpha,'bca',true,false);


        %% get column spacing
        getCI_PwCS =  false;
        [average_spacing_mm,local_spacing_mm,newROI] = getColumnsSpacing(data_obj,DataFolder,0.1,2,0.1,getCI_PwCS,true);

        %% get Pw infos
        lowpass_cutoffs_mm = 0.1:0.01:0.5;
        SizeGaussKernelPwDensityCalc = 0.5;
        do_plotting = false;
        PwInfo = getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI_PwCS,do_plotting,lowpass_cutoffs_mm,SizeGaussKernelPwDensityCalc);


%         %% find pws
%         [PwInfo.NumberPw,PwInfo.aniso,PwInfo.x_angle,PwInfo.PWxList,PwInfo.PWyList,PwInfo.signList, PwInfo.contours] = find_pinwheels(z,0,data_obj.ROI);
%         

        %% make map borders ROI
        ROI =data_obj.ROI;
        [YROI,XROI] = find(ROI);
        % [Xmin, Xmax] = findBorders(XROI);
        % [Ymin, Ymax] = findBorders(YROI);

        [x_range,y_range] = getRangeXY_ROI(data_obj_dir.ROI);

        average_spacing_px = average_spacing_mm*data_info.pix_per_mm; 

        %% get new x/y_range at the center of the ROI
        x_center = round(x_range(1) + (x_range(2)-x_range(1))/2);
        y_center = round(y_range(1) + (y_range(2)-y_range(1))/2);
        width_y_range = 3*average_spacing_px;
        width_x_range = 2*average_spacing_px;
        
        x_range = [x_center-width_x_range/2 x_center+width_x_range/2];
        y_range = [y_center-width_y_range/2 y_center+width_y_range/2];

        Xmin = x_range(1);
        Ymin = y_range(1);

        %% Plot Bloodvessels
        PlotBloodVessels(BloodVesselImg,ones(size(BloodVesselImg)),1)
        hold on
        contour(ROI,[1 1],'red','linewidth',2)
        hold on
        mm = 0.2;
        width_scale_pix = 5;
        spacing_pix = mm * data_info.pix_per_mm;
        plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-red')
        hold on
        text(Xmin+width_scale_pix,Ymin+width_scale_pix+2,[num2str(mm) ' mm'],'Color','red')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        xlim(x_range)
        ylim(y_range)

        hold on; set(gca,'view',[rotate rotate])
        print('-depsc', [DataFolder 'BloodVessels_' animal num2str(experiment_num) '.eps'])
        close

        %% plot polar map
        f1 = figure();
        plot_map(z,data_obj.ROI,3*sqrt(mean(abs(z(ROI)).^2)),1)
        xlim(x_range)
        ylim(y_range)
        title('orientation pref. map ')
        hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'PolarOri.ps'])

        %% plot ori map big
        f1=figure;
        plot_map(z,data_obj_dir.ROI,0,1)

        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','MarkerSize',5)
        hold on; SizesCI = getConfidenceRegionPw(PwInfo.pinwheel_stats,data_info.field_size_pix,1-alpha,true,false,false,'black',2);
        %getConfidenceRegionPw(pinwheel_stats,field_size_pix,Confidence,plot_CI,Plot_CI_Size_Scale,plotProb,color,linewidth)

        hold on; contour(NoiseMeasure.CI_ori,[30 30],'white','linewidth',2,'LineStyle','-')
        %hold on; contour(NoiseMeasure.CI_ori,[40 40],'white','linewidth',2,'LineStyle','--')

        xlim(x_range)
        ylim(y_range)
        title('orientation pref. map ')
        hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'Ori.ps'])

        %% plot ori map Abs big
        f1=figure;
        plot_Abs(abs(z).*ROI,gca)
        xlim(x_range)
        ylim(y_range)
        title('selectivity orientation map')
        hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'OriAbs.ps'])


        %% plot Pw CI
        f1 = figure;
        sigma = 0.05;
        PWxList = PwInfo.pinwheel_stats.x(~isnan(PwInfo.pinwheel_stats.x));
        PWyList = PwInfo.pinwheel_stats.y(~isnan(PwInfo.pinwheel_stats.y));
    
        local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PWxList, PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
        local_pw_dens = local_pw_dens./sum(local_pw_dens(ROI)).*sum(PwInfo.pinwheel_stats.probability,'all')/(sum(ROI,'all')/data_info.pix_per_mm^2);
    
        plot_mapAbs(local_pw_dens,'Pinwheel Prob. Density [1/mm^2]',max(local_pw_dens(data_obj.ROI),[],'all'),min(local_pw_dens(data_obj.ROI),[],'all'),data_obj.ROI,gca)
        hold on; contour(ROI,[1 1],'white','linewidth',2)
        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','linewidth',linewidth)
%         hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
%         hold on
%         contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
%         hold on
%         contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    %         hold on
    %         for ii = 1:size(PwInfo.pinwheel_stats.x,1)
    %             text(PwInfo.pinwheel_stats.x(ii,1)+2,PwInfo.pinwheel_stats.y(ii,1),num2str(PwInfo.pinwheel_stats.probability(ii)),'Color',color_contur,'FontSize',5)
    %         end
        xlim(x_range)
        ylim(y_range)
    
        hold on; set(gca,'view',[rotate rotate])
        yticks([])
        xticks([])
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'PwCI.ps'])


        %% plot Ori CI
        f1 = figure;
        plot_mapAbs(NoiseMeasure.CI_ori,'CI ori prefrerence [°]',180,0,ROI,gca)
        xlim(x_range)
        ylim(y_range)
        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','linewidth',linewidth)
        hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'OriCI.ps'])
        

        %% plot dir map big
        f1=figure;
        plot_map(z_dir,data_obj_dir.ROI,0,1)
        [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
        U=cos(angle(z_dir));V=sin(angle(z_dir));
        U(~ROI) = 0;
        V(~ROI) = 0;
        d=10;
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim(x_range)
        ylim(y_range)
        title('direction map ')
        hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'Dir.ps'])

        %% plot dir map Abs big
        f1=figure;
        plot_Abs(abs(z_dir).*ROI,gca)
        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'.','Color','white','linewidth',linewidth)
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)

        hold on
        contour(real(z_dir),[0 0],'linewidth',linewidth,'Color','white')
        hold on
        contour(imag(z_dir),[0 0],'linewidth',linewidth,'Color','white')

        xlim(x_range)
        ylim(y_range)
        title('selectivity direction map')
        hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirAbs.ps'])

        %% plot Dir CI map
        f1 = figure;
        plot_mapAbs(NoiseMeasure.CI_dir,'CI dir prefrerence [°]',360,0,ROI,gca)
        xlim(x_range)
        ylim(y_range)
        hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirCI.ps'])

        %% compare dir selectivity histogramms with norm pdf
        f1 = figure;
        histogram(abs(z_dir(ROI)),'Normalization','pdf','NumBins',20,'DisplayName','all sel.')
        hold on
        sel_pw = abs(z_dir(round(PwInfo.PWxList),round(PwInfo.PWyList)));
        histogram(sel_pw,'Normalization','pdf','NumBins',20,'DisplayName','pinwheel sel.')
        legend()
        print(f1,'-dpsc','-fillpage', [DataFolder animal num2str(experiment_num) 'DirSelHist.ps'])

%         [dir_average_spacing_mm,dir_local_spacing_mm,dir_newROI,dir_WavletCoefficient] = get_column_spacingManuel(z_dir,data_obj.ROI,data_obj.info.pix_per_mm,0.1,2,0.1);

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


    end
    

    
    
    
    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return 
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

