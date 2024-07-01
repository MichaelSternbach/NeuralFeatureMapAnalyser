function CheckDifferenzeMaps(animal,experiment_Num,AnimalDataFolder,DataFolderMain,FigureFolder)

    
    %% parameter
    experiment_Num = checkFormatNum(experiment_Num);
    
    
    rectangle = false;
    
    
    %% figure folder
    mkdir(FigureFolder)
    
    FigureFile = [FigureFolder animal '_CheckDifferenzeMaps'];

    
    rm_cmd = ['rm -f ' FigureFile '.ps'];
    disp(rm_cmd)
    system(rm_cmd)
    
    %% contour plots
    grey = [0.6 0.6 0.6];
    color_contur = grey;%'grey';
    linewidth = .4;
    set(gca,'Fontsize',20)
    
    if length(experiment_Num) == 1
        experiment_Num = 1: experiment_Num;
    end
    first = true;
    for experiment_num = experiment_Num
        
        %% data folder
        DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
        
        %% animal
        [data_info,~,data_obj,~,BloodVesselImg] = getAnimalData(animal,experiment_num,1,false,AnimalDataFolder);
        [average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,DataFolder,false);
        z = data_obj.filter_map(data_obj.read_map());
        data_info.ID = replace(data_info.ID,'_',' ');

        %% prepare rectangel

        if rectangle ~= false

            Full_ROI = data_obj.ROI;
            RectangleROI = getRectangleROI(data_info.rectangle,data_obj.ROI);
            data_obj.set_ROI(RectangleROI)
        end

        %% make map borders ROI
        [YROI,XROI] = find(data_obj.ROI);
        [Xmin, Xmax] = findBorders(XROI);
        [Ymin, Ymax] = findBorders(YROI);

        %% make Figure
        close all
        f = figure();
        t = tiledlayout(3,2);
        s=1;
        f.Position = [100 100 594*s 841*s];
        title(t,data_info.ID)
        %f.Position = [100 100 1000 8000];

        
        if string(lower(animal)) == "dunnart"
            mm=0.2;
        else
            mm=1;
        end
        
        width_scale_pix= 15;
        spacing_pix = mm * data_obj.info.pix_per_mm;
        
        %% Brain Surface
%         nexttile;
%         PlotBloodVessels(BloodVesselImg,ones(size(BloodVesselImg)),1)
%         hold on
%         contour(data_obj.ROI,[1 1],'red')
%         if rectangle ~= false
%             contour(Full_ROI,[1 1],'red')
%         end
% 
%         
%         hold on
%         plot([width_scale_pix,width_scale_pix+spacing_pix],[width_scale_pix, width_scale_pix],'-red')
%         hold on
%         text(width_scale_pix+spacing_pix+2,width_scale_pix,[num2str(mm) ' mm'],'Color','red')
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         title('Bloodvessel Map')

        %% difference Maps
        ax = nexttile;
        stim_orient = data_obj.info.stim_order;
        i1 = find(stim_orient==0);
        i2 = find(stim_orient==90);
        DiffernceMap = mean(data_obj.data(:,:,i1,:),4)-mean(data_obj.data(:,:,i2,:),4);
        if rectangle ~= false
            %plot_mapAbs(DiffernceMap,'Difference Map Cardinal',max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),Full_ROI,ax)
            PlotBloodVessels(DiffernceMap,Full_ROI,1)
            title('Difference Map Cardinal')
            contour(data_obj.ROI,[1 1],'red')
        else
            %plot_mapAbs(DiffernceMap,['Difference Map Cardinal'],max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),data_obj.ROI,ax)
            PlotBloodVessels(DiffernceMap,data_obj.ROI,1)
            title('Difference Map Cardinal')
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])

        hold on
        plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-red')
        hold on
        text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,[num2str(mm) ' mm'],'Color','red')
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])

        ax = nexttile;
        stim_orient = data_obj.info.stim_order;
        i1 = find(stim_orient==45);
        i2 = find(stim_orient==135);
        DiffernceMap = mean(data_obj.data(:,:,i1,:),4)-mean(data_obj.data(:,:,i2,:),4);
        if rectangle ~= false
            %plot_mapAbs(DiffernceMap,'Difference Map Oblique',max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),Full_ROI,ax)
            PlotBloodVessels(DiffernceMap,Full_ROI,1)
            title('Difference Map Oblique')
            contour(data_obj.ROI,[1 1],'red')
        else
            %plot_mapAbs(DiffernceMap,['Difference Map Oblique'],max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),data_obj.ROI,ax)
             PlotBloodVessels(DiffernceMap,data_obj.ROI,1)
            title('Difference Map Oblique')
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        
        
        %% difference Maps filtered
        ax = nexttile;
        stim_orient = data_obj.info.stim_order;
        i1 = find(stim_orient==0);
        i2 = find(stim_orient==90);
        DiffernceMap = mean(data_obj.data(:,:,i1,:),4)-mean(data_obj.data(:,:,i2,:),4);
        DiffernceMap = data_obj.filter_map(DiffernceMap);
        if rectangle ~= false
            %plot_mapAbs(DiffernceMap,'Difference Map Cardinal',max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),Full_ROI,ax)
            PlotBloodVessels(DiffernceMap,Full_ROI,1)
            title('Difference Map filtered')
            contour(data_obj.ROI,[1 1],'red')
        else
            %plot_mapAbs(DiffernceMap,['Difference Map Cardinal'],max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),data_obj.ROI,ax)
            PlotBloodVessels(DiffernceMap,data_obj.ROI,1)
            title('Difference Map filtered')
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])

        hold on
        plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-red')
        hold on
        text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,[num2str(mm) ' mm'],'Color','red')
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])

        ax = nexttile;
        stim_orient = data_obj.info.stim_order;
        i1 = find(stim_orient==45);
        i2 = find(stim_orient==135);
        DiffernceMap = mean(data_obj.data(:,:,i1,:),4)-mean(data_obj.data(:,:,i2,:),4);
        DiffernceMap = data_obj.filter_map(DiffernceMap);
        if rectangle ~= false
            %plot_mapAbs(DiffernceMap,'Difference Map Oblique',max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),Full_ROI,ax)
            PlotBloodVessels(DiffernceMap,Full_ROI,1)
            title('Difference Map filtered')
            contour(data_obj.ROI,[1 1],'red')
        else
            %plot_mapAbs(DiffernceMap,['Difference Map Oblique'],max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),data_obj.ROI,ax)
             PlotBloodVessels(DiffernceMap,data_obj.ROI,1)
            title('Difference Map filtered')
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        
        
        %% real and imaginary parts
        ax = nexttile;
        
        %plot_mapAbs(DiffernceMap,'Difference Map Cardinal',max(DiffernceMap(data_obj.ROI),[],'all'),min(DiffernceMap(data_obj.ROI),[],'all'),Full_ROI,ax)
        PlotBloodVessels(real(z),data_obj.ROI,1)
        title('Re OPM')

        set(gca,'xtick',[])
        set(gca,'ytick',[])

        hold on
        plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-red')
        hold on
        text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,[num2str(mm) ' mm'],'Color','red')
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])

        ax = nexttile;
        PlotBloodVessels(imag(z),data_obj.ROI,1)
        title('Im OPM')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])

        %% load pinwheel data
% 
%         getCI = false;
%         do_plotting=0;
%         llp_cutoffs = linspace(0.01, 1,100);
%         beta=0.5;

%         Bootstrapsamples = 100;
%         data_obj.prepare_samples_array(Bootstrapsamples)
        %PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta);
%         tracker_obj = pinwheel_tracker;
%         simple_track=true;
%         [PwInfo.pinwheel_stats,PwInfo.pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);


        
        %[~,~,~,PWxList,PWyList,~, ~] = find_pinwheels(z,0,data_obj.ROI);
        
%         %% plot OPMap
%         nexttile;
%         scale = (data_info.pix_per_mm*average_spacing_mm).^-1;
%         spacing_pix = 1/scale;
%         plot_map(z,data_obj.ROI,0,1)
%         hold on
%         plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-white')
%         hold on
%         text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,'Λ','Color','white')
%         hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color','white','linewidth',linewidth)
%         xlim([Xmin Xmax])
%         ylim([Ymin Ymax])


        
        %% save figure
        %print(f,'-depsc', [FigureFile '.eps'])
        %savefig(f,[FigureFile '.fig'])
        print(f, '-dpsc','-fillpage', '-append', [FigureFile '.ps'])
    end
    

    
    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return 
end







