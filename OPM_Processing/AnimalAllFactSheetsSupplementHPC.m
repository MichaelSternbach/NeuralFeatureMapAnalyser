function AnimalAllFactSheetsSupplementHPC(animal,experiment_Num,AnimalDataFolder,DataFolderMain,FigureFolder)

    
    %% parameter
    experiment_Num = checkFormatNum(experiment_Num);
    
    
    rectangle = false;
    
    
    %% figure folder
    mkdir(FigureFolder)
    
    if rectangle == false
         FigureFile = [FigureFolder animal '_FactSheet'];
    else
        FigureFile = [FigureFolder animal '_FactSheet_rectangle'];
    end
    
    rm_cmd = ['rm -f ' FigureFile '.ps'];
    disp(rm_cmd)
    system(rm_cmd)
    
    %% contour plots
    grey = [0.6 0.6 0.6];
    color_contur = grey;%'grey';
    linewidth = 2;
    
    
    for experiment_num = 1: experiment_Num
        
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

        f = figure();
        t = tiledlayout(5,3);
        title(t,data_info.ID)
        %f.Position = [100 100 1000 8000];

        %% Brain Surface
        if string(lower(animal)) == "dunnart"
            mm=0.2;
        else
            mm=1;
        end
        
        width_scale_pix= 15;

        nexttile;
        PlotBloodVessels(BloodVesselImg,ones(size(BloodVesselImg)),1)
        hold on
        contour(data_obj.ROI,[1 1],'red')
        if rectangle ~= false
            contour(Full_ROI,[1 1],'red')
        end

        spacing_pix = mm * data_obj.info.pix_per_mm;
        hold on
        plot([width_scale_pix,width_scale_pix+spacing_pix],[width_scale_pix, width_scale_pix],'-red')
        hold on
        text(width_scale_pix+spacing_pix+2,width_scale_pix,[num2str(mm) ' mm'],'Color','red')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        title('Bloodvessel Map')

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

        %% load pinwheel data

        getCI = false;
        do_plotting=0;
        llp_cutoffs = linspace(0.01, 1,100);
        beta=0.5;

        Bootstrapsamples = 100;
        data_obj.prepare_samples_array(Bootstrapsamples)
        PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta);
%         tracker_obj = pinwheel_tracker;
%         simple_track=true;
%         [PwInfo.pinwheel_stats,PwInfo.pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);


        
        %[~,~,~,PWxList,PWyList,~, ~] = find_pinwheels(z,0,data_obj.ROI);
        
        %% plot OPMap
        
        scale = (data_info.pix_per_mm*average_spacing_mm).^-1;
        spacing_pix = 1/scale;
        plot_map(z,data_obj.ROI,0,1)
        hold on
        plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-white')
        hold on
        text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,'Λ','Color','white')
        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color','white','linewidth',linewidth)
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])

        %% Wavelength map
        %[~,local_spacing_mm,~,~,~,~] =  loadColumnsSpacing(data_obj,DataFolder,true,true);

        ax = nexttile;
        plot_mapAbs(local_spacing_mm,['Local Column Spacing [mm]'],max(local_spacing_mm(data_obj.ROI),[],'all'),min(local_spacing_mm(data_obj.ROI),[],'all'),data_obj.ROI,ax)
        if rectangle ~= false
            xlim([data_info.rectangle(1) data_info.rectangle(3)])
            ylim([data_info.rectangle(2) data_info.rectangle(4)])
        else
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])

        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        
        
        %% Wavelength CI map
        CISpacingFile = [DataFolder 'CI_MapSpacing_' data_obj.info.ID '.mat'];
        load(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm')
        
        AbsCI_local_spacing_mm = abs(CI_local_spacing_mm(:,:,2)-CI_local_spacing_mm(:,:,1));

        ax = nexttile;
        plot_mapAbs(AbsCI_local_spacing_mm,['CI Local Column Spacing [mm]'],max(AbsCI_local_spacing_mm(data_obj.ROI),[],'all'),min(AbsCI_local_spacing_mm(data_obj.ROI),[],'all'),data_obj.ROI,ax)
        if rectangle ~= false
            xlim([data_info.rectangle(1) data_info.rectangle(3)])
            ylim([data_info.rectangle(2) data_info.rectangle(4)])
        else
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])

        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        
        %% Filter data
        DataAndFilterFile = [DataFolder 'FilterFile.mat'];

        if isfile(DataAndFilterFile) %&& false 
            load(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumber','power_profile')
        else
            lowpass_cutoffs = linspace(0.2*average_spacing_mm,average_spacing_mm*.8,200);
            profile_range_mm = [0.01 2]; 
            profile_step_mm = 0.01;
            
            power_profile = define_filter_settings(data_info,data_obj.ROI,data_obj.data,profile_range_mm,profile_step_mm);
            lowpass_cutoffs = lowpass_cutoffs(lowpass_cutoffs<(data_obj.filter_parameters.highpass*.9));
            filtersPwNumber = find_lowpassPwNumber(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
            save(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumber','power_profile')
        end

        %% CI Maps Angle
        CI = calcCIs(data_obj,0.05,true,DataFolder);
        ROI = data_obj.ROI;

        maxMap = 180;
        minMap = 0;
        plotAbsTile(CI.SE.CI_angle, 'CI angle SE [°]',maxMap,minMap,ROI)
        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        if rectangle ~= false
            xlim([data_info.rectangle(1) data_info.rectangle(3)])
            ylim([data_info.rectangle(2) data_info.rectangle(4)])
        else
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
        end

        plotAbsTile(CI.BCA.CI_angle, 'CI angle BCA [°]',maxMap,minMap,ROI)
        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        if rectangle ~= false
            xlim([data_info.rectangle(1) data_info.rectangle(3)])
            ylim([data_info.rectangle(2) data_info.rectangle(4)])
        else
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
        end
        
        %% powerspectrum OPM
        nexttile;

        plot(average_spacing_mm./power_profile.scale_mm,power_profile.values,'DisplayName','Power Profile Unfiltered Map');
        hold on
        plot([average_spacing_mm./data_obj.filter_parameters.lowpass average_spacing_mm./data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','Lowpass Cutoff')
        hold on
        plot([average_spacing_mm./data_obj.filter_parameters.highpass average_spacing_mm./data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','Highpass Cutoff')

        xlabel('Wavevector (1/Λ)')
        ylabel('Power')
        xlim([average_spacing_mm/2 average_spacing_mm./(data_obj.filter_parameters.lowpass*0.8)])
        ylim([min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
        %set(gca,'fontsize',15)
        axis square
        %axis tight
        legend('Location','northoutside')
        

        %% CI Maps Abs
        [preMax,OrderMax] = getOrder(max([CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
        [preMin,OrderMin] = getOrder(max(-[CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
        maxMap = ceil(preMax)*10^OrderMax;
        minMap = -ceil(preMin)*10^OrderMin;
        plotAbsTile(CI.SE.CI_Abs, 'CI abs SE',maxMap,minMap,ROI)
        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        if rectangle ~= false
            xlim([data_info.rectangle(1) data_info.rectangle(3)])
            ylim([data_info.rectangle(2) data_info.rectangle(4)])
        else
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
        end
        plotAbsTile(CI.BCA.CI_Abs, 'CI abs BCA',maxMap,minMap,ROI)
        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        if rectangle ~= false
            xlim([data_info.rectangle(1) data_info.rectangle(3)])
            ylim([data_info.rectangle(2) data_info.rectangle(4)])
        else
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
        end
        

        %% pw Plateau
        nexttile;
        plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2),'DisplayName','Pinwheel Data')
        hold on
        minPw = min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2))*0.9;
        maxPw = max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2))*1.1;
        plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[minPw maxPw],'DisplayName','Lowpass Cutoff')
        %min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all')
        %xlim([0.2 .8].*average_spacing_mm)
        ylim([minPw maxPw])
        xlabel('Scale [mm]')
        ylabel('# Pinwheels')
        axis square;

        %% pinwheel position variability and pw prob.
        %nexttile;
%         plot_map(data_obj.filter_map(data_obj.read_map()),data_obj.ROI,0,1);
%         hold
%         plotPinwheelStats(PwInfo.pinwheel_stats,data_obj.info.field_size_pix)
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         if rectangle ~= false
%             xlim([data_info.rectangle(1) data_info.rectangle(3)])
%             ylim([data_info.rectangle(2) data_info.rectangle(4)])
%         else
%             xlim([Xmin Xmax])
%             ylim([Ymin Ymax])
%         end

        sigma = 0.05;
        PWxList = PwInfo.pinwheel_stats.x(~isnan(PwInfo.pinwheel_stats.x));
        PWyList = PwInfo.pinwheel_stats.y(~isnan(PwInfo.pinwheel_stats.y));

        local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PWxList, PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
        local_pw_dens = local_pw_dens./sum(local_pw_dens,'all').*sum(PwInfo.pinwheel_stats.probability,'all');

        ax = nexttile;
        plot_mapAbs(local_pw_dens,'Pinwheel Prob. Density',max(local_pw_dens(data_obj.ROI),[],'all'),min(local_pw_dens(data_obj.ROI),[],'all'),data_obj.ROI,ax)
        
        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        for ii = 1:size(PwInfo.pinwheel_stats.x,1)
            text(PwInfo.pinwheel_stats.x(ii,1)+2,PwInfo.pinwheel_stats.y(ii,1),num2str(PwInfo.pinwheel_stats.probability(ii)),'Color',color_contur,'FontSize',18)
        end

        xlim([Xmin Xmax])
        ylim([Ymin Ymax])

        %hold on; set(gca,'view',[rotate rotate])
        yticks([])
        xticks([])

        %% load pinwheel CI data
        CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
        load(CIPwFile,'CI_PwDensities','alpha','PwInfosBS','PwInfosJS')
        LocalPwDensity = CI_PwDensities.LocalPwDensityPlateuFit(:,:,2);
        CI_LocalPwDensity = abs(CI_PwDensities.LocalPwDensityPlateuFit(:,:,3)-CI_PwDensities.LocalPwDensityPlateuFit(:,:,1));
        
        %% plot local pw density
        ax = nexttile;
        plot_mapAbs(LocalPwDensity,'Local Pinwheel Density [1/Λ]',max(LocalPwDensity(data_obj.ROI),[],'all'),min(LocalPwDensity(data_obj.ROI),[],'all'),data_obj.ROI,ax)
        if rectangle ~= false
            xlim([data_info.rectangle(1) data_info.rectangle(3)])
            ylim([data_info.rectangle(2) data_info.rectangle(4)])
        else
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])

        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        
        %% plot CI local pw density
        ax = nexttile;
        plot_mapAbs(CI_LocalPwDensity,'CI Local Pinwheel Density [1/Λ]',max(CI_LocalPwDensity(data_obj.ROI),[],'all'),min(CI_LocalPwDensity(data_obj.ROI),[],'all'),data_obj.ROI,ax)
        if rectangle ~= false
            xlim([data_info.rectangle(1) data_info.rectangle(3)])
            ylim([data_info.rectangle(2) data_info.rectangle(4)])
        else
            xlim([Xmin Xmax])
            ylim([Ymin Ymax])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])

        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        
        %% make animal table
        Animal = {data_info.ID};
        MeanSpacing_mm = {num2str([CI_average_spacing_mm(1) average_spacing_mm CI_average_spacing_mm(2)])};
%         MeanPwDensity = {num2str(PwInfo.MeanPwDensity)};
        MeanPwDensityManuel = {num2str([CI_PwDensities.PwDensityPosEstimate(1) PwInfo.PwDensityPosEstimate CI_PwDensities.PwDensityPosEstimate(2)])};
        PlateauFitPwDensity = {num2str([CI_PwDensities.PwDensityPlateuFit(1) PwInfo.PwDensityPlateuFit CI_PwDensities.PwDensityPlateuFit(2)])};
        
        
        MeanPwDensityCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'MeanPwDensity',false,alpha);
        MeanPwDensity = {num2str([MeanPwDensityCI(1) PwInfo.MeanPwDensity MeanPwDensityCI(2)])};
        
        NumberPwCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'NumberPw',false,alpha);
        NumberPw = {num2str([NumberPwCI(1) PwInfo.NumberPw NumberPwCI(2)])};
        
        T = table(MeanSpacing_mm,NumberPw,MeanPwDensity,MeanPwDensityManuel,PlateauFitPwDensity,'RowNames',Animal);
        
        if experiment_num ==1
            AllAnimalTable = T;
        else
            AllAnimalTable = [AllAnimalTable; T];
        end
        
        T = rows2vars(T);
        T.Properties.VariableNames(1) = "Name";
        
%         %% plot single animal table
%         ax = nexttile;
%         uitable(ax,'Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
%         'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
%         
%         
        
        
        %% save figure
        %print(f,'-depsc', [FigureFile '.eps'])
        %savefig(f,[FigureFile '.fig'])
        print(f, '-dpsc','-fillpage', '-append', [FigureFile '.ps'])
    end
    
    %% print all animal table to seperate page
    disp(AllAnimalTable)
    
    fT= figure();
    
    uitable('Data',AllAnimalTable{:,:},'ColumnName',AllAnimalTable.Properties.VariableNames,...
    'RowName',AllAnimalTable.Properties.RowNames,'Units', 'Normalized'); %, 'Position',[0, 0, 1, 1]
    

    print(fT, '-dpsc','-fillpage', '-append', [FigureFile '.ps'])
    
    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return 
end







