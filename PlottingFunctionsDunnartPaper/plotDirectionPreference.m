function plotDirectionPreference(experiment_num_list,animal,AnimalDataFolder,FigureFolder,...
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
    linewidth = 1;
    Bootstrapsamples = 100;
%     
%     %% set filter settings
%     lowpass_mm = 0.25;
%     highpass_mm = 0.56;
    
    %% figure folder and files
    if ~isfolder(FigureFolder)
        mkdir(FigureFolder) 
    end

    FigureFile = [FigureFolder animal '_DirectionPrefernceGIF' num2str(GIF_SN_TH)];   
    
    %% remove old figure file
    if isfile([FigureFile '.ps'])
        delete([FigureFile '.ps'])
    end
  
%     rotate = 0;% -90;
    for experiment_num = experiment_num_list
        
        %% get data infos
%         data_info_file_list = '/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/OPM_Processing/experiment_info.csv';
%         [data_info,data_path] = info_handle(animal,experiment_num,AnimalDataFolder,data_info_file_list);
%         if isfield(data_info,'pix_per_mm')
%             data_info.pixels_per_mm = data_info.pix_per_mm;
%         else
%             data_info.pix_per_mm = data_info.pixels_per_mm;
%         end
        
        %% get data
        [data_info,data_path,data_obj,data,~] = getAnimalData(animal,experiment_num,AnimalDataFolder,data_info_file_list);
%         getDirectionData = false;
%         NonStimuliResponseRemoval = 'average';
%         data = NoAveragePreProcessRawDataJason(data_info.expIds,data_info.refWin,data_info.sigWin,data_info.partId,data_path,data_info.ID,getDirectionData,NonStimuliResponseRemoval);
%         
%         %% get ROI
%         info_path = [data_path,'exp_info.mat'];
%         load(info_path,'ROI');
%         
%         %% prepare data obj
%         data_obj = data_handle_corrected(data_info,data,ROI);
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

        %% make map borders ROI
        ROI =data_obj.ROI;
        [YROI,XROI] = find(ROI);
        [Xmin, Xmax] = findBorders(XROI);
        [Ymin, Ymax] = findBorders(YROI);
        
        %% get OPM map
        z = data_obj.filter_map(data_obj.read_map());
%         %data = data - mean(data(:,:,1:end-1,:),3:4);
%         data = data_obj.data;
%         StimDir = data_info.stim_order;
%         StimDir = StimDir(~isnan(StimDir));
%         StimDir = reshape(StimDir,[1 1 length(StimDir) 1]);
%         if GIF_SN_TH>0
%             data(:,:,1:end-1,:) = GIF(data(:,:,1:end-1,:),GIF_SN_TH);
%         end
%         z = mean(data(:,:,1:length(StimDir),:).*exp(2i*pi*StimDir/360),3:4);
%         z = data_obj.filter_map(z);

        %% find pws
        [PwInfo.NumberPw,PwInfo.aniso,PwInfo.x_angle,PwInfo.PWxList,PwInfo.PWyList,PwInfo.signList, PwInfo.contours] = find_pinwheels(z,0,data_obj.ROI);
        

        %% get dir map
%         getDirectionData = true;
%         DirectionData = NoAveragePreProcessRawDataJason(data_info.expIds,data_info.refWin,data_info.sigWin,data_info.partId,data_path,data_info.ID,getDirectionData,NonStimuliResponseRemoval);
        
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

        
        %% plot dir map big
        f1=figure;
        z_dir = data_obj_dir.filter_map(data_obj_dir.read_map(1,false,true));
        plot_map(z_dir,data_obj_dir.ROI,0,1)
        [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
        U=cos(angle(z_dir));V=sin(angle(z_dir));
        U(~ROI) = 0;
        V(~ROI) = 0;
        d=10;
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('direction map ')
%         hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [FigureFolder animal num2str(experiment_num) 'Dir.ps'])

        %% plot dir map Abs big
        f1=figure;
        plot_Abs(abs(z_dir).*ROI,gca)
        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('selectivity direction map')
%         hold on; set(gca,'view',[rotate rotate])
        print(f1,'-dpsc','-fillpage', [FigureFolder animal num2str(experiment_num) 'DirAbs.ps'])


        %% make Figure
        f = figure(1);
        t = tiledlayout(3,3);
        title(t,data_info.ID)
            
        
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


        %% plot dir map data obj
        z_dir = data_obj_dir.filter_map(data_obj_dir.read_map(1,false,true));
        nexttile(t)
        plot_map(z_dir,data_obj_dir.ROI,0,1)
        [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
        U=cos(angle(z_dir));V=sin(angle(z_dir));
        U(~ROI) = 0;
        V(~ROI) = 0;

        d=10;
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('direction map ')

        %% plot ori map
        nexttile(t)
        plot_map(z,data_obj.ROI,0,1)
        [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
        U=cos(angle(z_dir));V=sin(angle(z_dir));
        U(~ROI) = 0;
        V(~ROI) = 0;

        d=10;
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('orientation map')




        %% plot ori map data obj
        z_av = data_obj_dir.filter_map(data_obj_dir.read_map());
        nexttile(t)
        plot_map(z_av,data_obj.ROI,0,1)
        [X, Y]= meshgrid(1:size(z_dir,2),1:size(z_dir,1));
        U=cos(angle(z_dir));V=sin(angle(z_dir));
        U(~ROI) = 0;
        V(~ROI) = 0;

        d=10;
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('orientation map av dir')

        [PwInfo_av.NumberPw,PwInfo_av.aniso,PwInfo_av.x_angle,PwInfo_av.PWxList,PwInfo_av.PWyList,PwInfo_av.signList, PwInfo_av.contours] = find_pinwheels(z_av,0,data_obj_dir.ROI);


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
% 
% 
%         %% plot CI direction map
%         [x_range,y_range] = getRangeXY_ROI(data_obj_dir.ROI);
%         alpha = 0.05;
%         [NoiseMeasure.CI_dir,~,~] = getCI(data_obj_dir,alpha,'bca',true,true);
%         ax = nexttile(t);
%         plot_mapAbs(NoiseMeasure.CI_dir,'CI dir prefrerence [°]',360,0,ROI,ax)
%         xlim(x_range)
%         ylim(y_range)
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



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               
        %% plot polar map
        nexttile(t)
        plot_map(z_dir,ROI,3*sqrt(mean(abs(z_dir(ROI)).^2)),1)
        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('polar direction map')

        %% plot selectivity direction map
        ax = nexttile(t);
        plot_Abs(abs(z_dir).*ROI,ax)
        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('selectivity direction map')

        %% plot selectivity orientation map
        ax = nexttile(t);
        plot_Abs(abs(z).*ROI,ax)
        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        % hold on
        % quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('selectivity orientation map')
    


        %% plot polar direction map qith average over opposite directions pinwheel centers
        nexttile(t)
        plot_map(z_dir,ROI,3*sqrt(mean(abs(z_dir(ROI)).^2)),1)
        hold on; plot(PwInfo_av.PWxList,PwInfo_av.PWyList,'x','Color','white','linewidth',linewidth)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('polar direction map av Pws')

        %% plot selectivity direction map with average over opposite directions pinwheel centers
        ax = nexttile(t);
        plot_Abs(abs(z_dir).*ROI,ax)
        hold on; plot(PwInfo_av.PWxList,PwInfo_av.PWyList,'x','Color','white','linewidth',linewidth)
        hold on
        quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('selectivity direction map av Pws')

        %% plot selectivity orientation map with average over opposite directions pinwheel centers
        ax = nexttile(t);
        plot_Abs(abs(z_av).*ROI,ax)
        hold on; plot(PwInfo_av.PWxList,PwInfo_av.PWyList,'x','Color','white','linewidth',linewidth)
        % hold on
        % quiver(X(1:d:end,1:d:end),Y(1:d:end,1:d:end),U(1:d:end,1:d:end),V(1:d:end,1:d:end),0.3,'-white', 'LineWidth', 0.9)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('selectivity orientation map av ')

        % %% get response data 
        % Filter = true;
        % ResponseDataDir = getResponseData(DirectionData,Filter,data_obj);
        % ResponseDataOri = getResponseData(data,Filter,data_obj);

        % %% plot responses

        % direction = 5;
        % period = 8;
        % ax=nexttile(t);
        % plot_Abs(ResponseDataOri(:,:,direction).*ROI,ax)
        % hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        % xlim([Xmin Xmax])
        % ylim([Ymin Ymax])
        % title(['orientation' num2str(data_info.stim_order(direction)) '°'])

        % ax = nexttile(t);
        % plot_Abs(mean(ResponseDataDir(:,:,[direction direction+period]),3).*ROI,ax)
        % hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        % xlim([Xmin Xmax])
        % ylim([Ymin Ymax]) 
        % title(['mean direction' num2str(data_info.stim_order(direction)) '° and ' num2str(data_info.stim_order(direction)+180) '°'])

        % ax=nexttile(t);
        % plot_Abs(ResponseDataDir(:,:,direction,:).*ROI,ax)
        % hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        % xlim([Xmin Xmax])
        % ylim([Ymin Ymax])
        % title(['direction' num2str(data_info.stim_order(direction))  '°'])
       
        % ax=nexttile(t);
        % plot_Abs(ResponseDataDir(:,:,direction+period).*ROI,ax)
        % hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        % xlim([Xmin Xmax])
        % ylim([Ymin Ymax])
        % title(['direction'  num2str(data_info.stim_order(direction)+180) '°'])

        % ax=nexttile(t);
        % plot_Abs(ResponseDataDir(:,:,direction)-ResponseDataDir(:,:,direction+period).*ROI,ax)
        % hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'x','Color','white','linewidth',linewidth)
        % xlim([Xmin Xmax])
        % ylim([Ymin Ymax])

        % title(['difference direction' num2str(data_info.stim_order(direction)) '° and ' num2str(data_info.stim_order(direction)+180) '°'])
        
        print(f, '-dpsc','-fillpage','-append', [FigureFile '.ps'])%'-append',
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

