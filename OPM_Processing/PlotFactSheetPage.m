function [AnimalTable,ResultData] = PlotFactSheetPage(data_info,data_obj,DataFolder,FigureFile,...
    BiasDataFolder,width_scale_pix,color_contur,linewidth,Fontsize,alpha)


    
    %% set default values
    if nargin < 5
        FigureFile = [DataFolder 'Factsheet_' data_info.animal ' ' data_info.ID '.ps'];
    end
    if nargin < 6
        BiasDataFolder = '';
        %%BiasDataFolder = ['~/Cloud/Cloud/PhD/Writing/phd_thesis/OPM_Methods/MakeNoiseFromDataROI_ColumnSpacing/ResultDatav0/' animal '/'];
    end
    if nargin < 7
        width_scale_pix = 15;
    end
    if nargin < 8
        color_contur = [0.6 0.6 0.6];%%'grey';
    end
    if nargin < 9
        linewidth = 0.3;
    end
    if nargin < 10
        Fontsize = 20;
    end
    set(gca,'Fontsize',Fontsize)
    if nargin <11
        alpha = 0.05;
    end

    if alpha ~= 0.3180
        CI_name = [num2str(round((1-alpha)*100,1)) '% CI '];
    else
        CI_name = 'SE ';
    end


    % if string(lower(animal)) == "dunnart"
    %     mm=0.2;
    % else
    %     mm=1;
    % end

    %% get map

    z = data_obj.filter_map(data_obj.read_map());
    
    
    %% get column spacingmap and ID
%     [data_info,~,data_obj,~,BloodVesselImg] = getAnimalData(animal,experiment_num,AnimalDataFolder);
    
    %% load spacing data
    SpacingFile = [DataFolder 'MapSpacingFiltered_' data_info.ID '.mat'];
    load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI')
    


    %% load pinwheel data

    getCI_ = false;
    do_plotting=0;
    llp_cutoffs = linspace(0.01, 1,100);
    beta=0.5;

    Bootstrapsamples = 100;
    data_obj.prepare_samples_array(Bootstrapsamples)
    PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI_,do_plotting,llp_cutoffs,beta);
    
%     %% load CI data pw
% 
%     CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
%     disp('load BS and JS')
%     load(CIPwFile,'PwInfosBS','PwInfosJS')
%     
%     disp('calc Pw CIs fro BS and JS')
%     CI_PwDensities.PwDensityPosEstimate = getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'PwDensityPosEstimate',false,alpha);
%     CI_PwDensities.PwDensityPlateuFit = getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'PwDensityPlateuFit',false,alpha);
%     CI_PwDensities.LocalPwDensityPlateuFit = getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'automated_pw_density',true,alpha);
%     
    %% load CI data spacing
    CISpacingFile = [DataFolder 'CI_MapSpacingFiltered_' data_obj.info.ID '.mat'];
    disp('load BS and JS')
    load(CISpacingFile,'average_spacings_mm','local_spacings_mm','jackstat_average_spacing_mm','local_spacingsJS_mm','newROIsBS','newROIsJS')
    
    disp('calc Pw CIs fro BS and JS')

    bootstat_local_spacings_mm = convert2vector(local_spacings_mm,data_obj,2);
    jackstat_local_spacing_mm =  convert2vector(local_spacingsJS_mm,data_obj,1);

    CI_average_spacing_mm = bootstrap_ci(average_spacings_mm(2:end),average_spacing_mm,jackstat_average_spacing_mm,alpha);
    CI_local_spacing_mmVector = bootstrap_ci(bootstat_local_spacings_mm(:,2:end),data_obj.array2vector(local_spacing_mm),jackstat_local_spacing_mm,alpha);           
    CI_local_spacing_mm = zeros([size(data_obj.ROI) 2]);
    CI_local_spacing_mm(:,:,1) = data_obj.vector2array(CI_local_spacing_mmVector(:,1));
    CI_local_spacing_mm(:,:,2) = data_obj.vector2array(CI_local_spacing_mmVector(:,2));   


%     %% load spacing data
%     SpacingFile = [DataFolder 'MapSpacingFiltered_' data_info.ID '.mat'];
%     load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI')
% 
%     %% load CI spacing data
%     CISpacingFile = [DataFolder 'CI_MapSpacingFiltered_' data_info.ID '.mat'];
%     %load(CISpacingFile,'jackstat_average_spacing_mm','average_spacings_mm')
%     load(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm','average_spacings_mm','local_spacings_mm','newROIs','jackstat_average_spacing_mm')
%     %CI_average_spacing_mm = bootstrap_ci(average_spacings_mm(2:end),average_spacing_mm,jackstat_average_spacing_mm,alpha);


    %% load pinwheel data
    PwInfoFile = [DataFolder 'PwInfo_' data_info.ID '.mat'];
    load(PwInfoFile,'PwInfo')


    %% load pinwheel CI data
    CIPwFile = [DataFolder 'CI_PwDensity_' data_info.ID '.mat'];
    load(CIPwFile,'alpha','PwInfosBS','PwInfosJS')
    %load(CIPwFile,'CI_PwDensities')

    %% get CI data ori
    DoFilter= true;
    [CI.BCA.CI_angle,CI.BCA.CI_Abs,CI.BCA.ROI] = getCI(data_obj,alpha,'bca',DoFilter);
    ROI = data_obj.ROI;

    %% adapt ID for plotting
    data_info.ID = replace(data_info.ID,'_',' ');

    %% make map borders ROI
    [YROI,XROI] = find(data_obj.ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);


    %% Borders Abs (CI) Plots
    [preMax,OrderMax] = getOrder(max([CI.BCA.CI_Abs ],[],'all'));%CI.SE.CI_Abs
    [preMin,OrderMin] = getOrder(max(-[CI.BCA.CI_Abs],[],'all'));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMin)*10^OrderMin;


    %% make Figure
    close all
    f = figure();
    t = tiledlayout(4,3);
    s=1;
    f.Position = [100 100 594*s 841*s];
    title(t,data_info.ID)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% plot OPMap
    
    nexttile;
    scale = (data_info.pix_per_mm*average_spacing_mm).^-1;
    spacing_pix = 1/scale;
    plot_map(z./abs(z),data_obj.ROI,0,1)
    hold on
    plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-white')
    hold on
    text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,'Λ','Color','white')
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color','white','linewidth',linewidth)
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
    title('Orientation Preference')

    

    
    %% CI Maps Angle
    maxMap = 180;
    minMap = 0;
    
    

    plotAbsTile(CI.BCA.CI_angle, [CI_name 'Prefernce [°]'],maxMap,minMap,ROI)
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
    hold on
    contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    hold on
    contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)

    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
    
    %% plot CPDF CI Maps Angle
    ax = nexttile;
    plotCPDF(CI.BCA.CI_angle(ROI),'','-',ax)
    title([CI_name 'Prefernce CPDF'])
    xlabel([CI_name 'Prefernce ≤ X [°]'])
    ylabel('% of pix.')
    %xlim([0,1])
    axis(ax,'square')
    ResultData.CI_angle = CI.BCA.CI_angle(ROI);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% plot polar map
    nexttile;
    scale = (data_info.pix_per_mm*average_spacing_mm).^-1;
    spacing_pix = 1/scale;
    plot_map(z,data_obj.ROI,3*sqrt(mean(abs(z(ROI)).^2)),1)
    hold on
    plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-white')
    hold on
    text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,'Λ','Color','white')
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color','white','linewidth',linewidth)
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
    title('Polar Map')

    %% plot selectivity Abs
    maxMap=max(abs(z(ROI)),[],'all');%7;
    plotAbsTile(abs(z), 'selectivity',maxMap,minMap,ROI)
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
    hold on
    contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    hold on
    contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])

    
%         %% CI Maps Abs
%         
%         plotAbsTile(CI.BCA.CI_Abs, 'CI abs BCA',maxMap,minMap,ROI)
%         hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
%         hold on
%         contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
%         hold on
%         contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
% 
%         xlim([Xmin Xmax])
%         ylim([Ymin Ymax])

    
    %% plot CPDF CI Maps Abs
    ax = nexttile;
    CIAbs= CI.BCA.CI_Abs./(2*abs(z));
    plotCPDF(CIAbs(ROI),'','-',ax)
    title(['rel. ' CI_name 'Selectivity CPDF'])
    xlabel([CI_name 'Selectivity ≤ X [Selectivity]'])
    ylabel('% of pix.')
    %xlim([0,1])
    axis(ax,'square')
    ResultData.CI_Abs = CI.BCA.CI_Abs(ROI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% pinwheel position variability 

    sigma = 0.05;
    PWxList = PwInfo.pinwheel_stats.x(~isnan(PwInfo.pinwheel_stats.x));
    PWyList = PwInfo.pinwheel_stats.y(~isnan(PwInfo.pinwheel_stats.y));

    local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PWxList, PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
    local_pw_dens = local_pw_dens./sum(local_pw_dens,'all').*sum(PwInfo.pinwheel_stats.probability,'all');

    ax = nexttile;
    plot_mapAbs(local_pw_dens,'Pinwheel Prob. Density',max(local_pw_dens(data_obj.ROI),[],'all'),min(local_pw_dens(data_obj.ROI),[],'all'),data_obj.ROI,ax)
    
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
    hold on
    contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    hold on
    contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
%         hold on
%         for ii = 1:size(PwInfo.pinwheel_stats.x,1)
%             text(PwInfo.pinwheel_stats.x(ii,1)+2,PwInfo.pinwheel_stats.y(ii,1),num2str(PwInfo.pinwheel_stats.probability(ii)),'Color',color_contur,'FontSize',5)
%         end

    xlim([Xmin Xmax])
    ylim([Ymin Ymax])

    %hold on; set(gca,'view',[rotate rotate])
    yticks([])
    xticks([])
    
    %% pw prob. and confidence region

%         sigma = 0.05;
%         PWxList = PwInfo.pinwheel_stats.x(~isnan(PwInfo.pinwheel_stats.x));
%         PWyList = PwInfo.pinwheel_stats.y(~isnan(PwInfo.pinwheel_stats.y));
% 
%         local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PWxList, PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
%         local_pw_dens = local_pw_dens./sum(local_pw_dens,'all').*sum(PwInfo.pinwheel_stats.probability,'all');
% 
    ax = nexttile;
    %plot_mapAbs(zeros(size(local_pw_dens)),'Pinwheel Prob. Density',max(local_pw_dens(data_obj.ROI),[],'all'),min(local_pw_dens(data_obj.ROI),[],'all'),data_obj.ROI,ax)
    imagesc(ax,zeros(size(ROI)))
    hold on;
    SizesCI = getConfidenceRegionPw(PwInfo.pinwheel_stats,data_info.field_size_pix,1-alpha);
    
%         hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
%         hold on
%         contour(real(z),[0 0],'linewidth',linewidth,'Color','white')%color_contur
%         hold on
%         contour(imag(z),[0 0],'linewidth',linewidth,'Color','white')
%         hold on
%         for ii = 1:size(PwInfo.pinwheel_stats.x,1)
%             text(PwInfo.pinwheel_stats.x(ii,1)+2,PwInfo.pinwheel_stats.y(ii,1),num2str(PwInfo.pinwheel_stats.probability(ii)),'Color',color_contur,'FontSize',5)
%         end

    contour(ROI,[1 1],'white','linewidth',linewidth)
    m=100;
    cm_viridis=viridis(m);
    colormap(ax, cm_viridis);
    axis image
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
    title([CI_name 'Pinwheel Positions'])
    %hold on; set(gca,'view',[rotate rotate])
    yticks([])
    xticks([])
    
%         %% plot CPDF pinwheel prob.
%         ax = nexttile;
%         plotCPDF(PwInfo.pinwheel_stats.probability,'','-',ax)
%         title('Pinwheel Prob.')
%         xlabel('PW prob. ≤ X')
%         ylabel('% of pinwheels')
%         xlim([0,1])
%         axis(ax,'square')
    
    
            
    %% plot CPDF pinwheel CI Size
    ax = nexttile;
    PwCI = sqrt(SizesCI/(data_info.pix_per_mm)^2/pi);
    plotCPDF(PwCI,'','-',ax)
    title(['Pinwheel ' CI_name 'Size CPDF'])
    xlabel(['PW ' CI_name 'eff. r ≤ x [mm]'])
    ylabel('% of pinwheels')
    %xlim([0,1])
    axis(ax,'square')
    ResultData.PwCI = PwCI;
    

%         %% pw Plateau
%         nexttile;
%         plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2),'DisplayName','Pinwheel Data')
%         hold on
%         minPw = min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2))*0.9;
%         maxPw = max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2))*1.1;
%         plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[minPw maxPw],'DisplayName','Lowpass Cutoff')
%         %min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all')
%         %xlim([0.2 .8].*average_spacing_mm)
%         ylim([minPw maxPw])
%         xlabel('Scale [mm]')
%         ylabel('# Pinwheels')
%         axis square;
% 
    %% load pinwheel CI data
    % LocalPwDensity = CI_PwDensities.LocalPwDensityPlateuFit(:,:,2);
    % CI_LocalPwDensity = abs(CI_PwDensities.LocalPwDensityPlateuFit(:,:,3)-CI_PwDensities.LocalPwDensityPlateuFit(:,:,1));
%         
%         %% plot local pw density
%         ax = nexttile;
%         plot_mapAbs(LocalPwDensity,'Local Pw Density [1/Λ²]',max(LocalPwDensity(data_obj.ROI),[],'all'),min(LocalPwDensity(data_obj.ROI),[],'all'),data_obj.ROI,ax)
%         if rectangle ~= false
%             xlim([data_info.rectangle(1) data_info.rectangle(3)])
%             ylim([data_info.rectangle(2) data_info.rectangle(4)])
%         else
%             xlim([Xmin Xmax])
%             ylim([Ymin Ymax])
%         end
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
% 
%         hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
%         hold on
%         contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
%         hold on
%         contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
%         
%         %% plot CI local pw density
%         ax = nexttile;
%         plot_mapAbs(CI_LocalPwDensity,'CI Local Pw Density [1/Λ²]',max(CI_LocalPwDensity(data_obj.ROI),[],'all'),min(CI_LocalPwDensity(data_obj.ROI),[],'all'),data_obj.ROI,ax)
%         if rectangle ~= false
%             xlim([data_info.rectangle(1) data_info.rectangle(3)])
%             ylim([data_info.rectangle(2) data_info.rectangle(4)])
%         else
%             xlim([Xmin Xmax])
%             ylim([Ymin Ymax])
%         end
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
% 
%         hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
%         hold on
%         contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
%         hold on
%         contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Wavelength map

    ax = nexttile;
    plot_mapAbs(local_spacing_mm,['Local Column Spacing [mm]'],0.6,0.2,data_obj.ROI,ax)
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])

    set(gca,'xtick',[])
    set(gca,'ytick',[])

    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
    hold on
    contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    hold on
    contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    
    
    %% Wavelength CI map   
    AbsCI_local_spacing_mm = abs(CI_local_spacing_mm(:,:,2)-CI_local_spacing_mm(:,:,1));

    ax = nexttile;
    %plot_mapAbs(AbsCI_local_spacing_mm,['CI Local Column Spacing [mm]'],max(AbsCI_local_spacing_mm(data_obj.ROI),[],'all'),min(AbsCI_local_spacing_mm(data_obj.ROI),[],'all'),data_obj.ROI,ax)
    plot_mapAbs(AbsCI_local_spacing_mm,[CI_name 'Local Column Spacing [mm]'],0.5,0,data_obj.ROI,ax)
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])

    set(gca,'xtick',[])
    set(gca,'ytick',[])

    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
    hold on
    contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    hold on
    contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)

    
    %% plot CPDF pinwheel CI Size
    ax = nexttile;
    CI_CS = AbsCI_local_spacing_mm./(local_spacing_mm*2);
    plotCPDF(CI_CS(ROI),'','-',ax)
    title(['rel. ' CI_name 'ColumnSpacing CPDF'])
    xlabel(['rel. ' CI_name '\Lambda/ ≤ X '])
    ylabel('% of pix')
    %xlim([0,1])
    axis(ax,'square')
    ResultData.CI_CS = CI_CS(ROI);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save figure
    %print(f,'-depsc', [FigureFile '.eps'])
    %savefig(f,[FigureFile '.fig'])
    print(f, '-dpsc','-fillpage', '-append', [FigureFile '.ps'])



    
    %% make animal table
    Animal = {data_info.ID};
    MeanSpacing_mm = {printValueCI(CI_average_spacing_mm(1),average_spacing_mm,CI_average_spacing_mm(2),2)};%{num2str([CI_average_spacing_mm(1) average_spacing_mm CI_average_spacing_mm(2)])};
%         MeanPwDensity = {num2str(PwInfo.MeanPwDensity)};
    %MeanPwDensityManuel = {num2str([CI_PwDensities.PwDensityPosEstimate(1) PwInfo.PwDensityPosEstimate CI_PwDensities.PwDensityPosEstimate(2)])};
    %PlateauFitPwDensity = {num2str([CI_PwDensities.PwDensityPlateuFit(1) PwInfo.PwDensityPlateuFit CI_PwDensities.PwDensityPlateuFit(2)])};
    NumHypercolumns = {PwInfo.NumHypercolumns};
    
    MeanPwDensityCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'WeightedPwDensityFixedFilter',false,alpha);
    MeanPwDensity = {printValueCI(MeanPwDensityCI(1), PwInfo.WeightedPwDensityFixedFilter, MeanPwDensityCI(2),2)};%{num2str([MeanPwDensityCI(1) PwInfo.MeanPwDensity MeanPwDensityCI(2)])};
    
    NumberPwCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'NumberPw',false,alpha);
    NumberPw = {printValueCI(NumberPwCI(1),PwInfo.NumberPw,NumberPwCI(2),0)};%{num2str([NumberPwCI(1) PwInfo.NumberPw NumberPwCI(2)])};
    
    RatioPrefCI = {round(sum(CI.BCA.CI_angle(ROI)<40,'all')/sum(ROI,'all')*100,1)};
    RatioSelectivityCI = {round(sum(CIAbs<0.5,'all')/size(CIAbs,1)/size(CIAbs,2)*100,1)};
    RatioPwProb = {round(sum(PwInfo.pinwheel_stats.probability>0.80,'all')/length(PwInfo.pinwheel_stats.probability)*100,1)};
    RatioPwCI ={round(sum(PwCI<0.1,'all')/size(PwInfo.pinwheel_stats.probability,1)/size(PwInfo.pinwheel_stats.probability,2)*100,1)};
    RatioColumnSpacingCI = {round(sum(CI_CS(ROI)<0.2,'all')/sum(ROI,'all')*100,1)};
    
    if ~isempty(BiasDataFolder)
        % load([BiasDataFolder 'PwDensityBiasCorrection.mat'],'PwDensityBiasCorrection')
        % MeanPwDensityCorrected = {printValueCI(MeanPwDensityCI(1)-PwDensityBiasCorrection(experiment_num), PwInfo.WeightedPwDensityFixedFilter-PwDensityBiasCorrection(experiment_num), MeanPwDensityCI(2)-PwDensityBiasCorrection(experiment_num),2)};
        % AnimalTable = table(MeanSpacing_mm,NumHypercolumns,NumberPw,MeanPwDensity,MeanPwDensityCorrected,RatioPrefCI,RatioSelectivityCI,RatioPwProb,RatioPwCI,RatioColumnSpacingCI,'RowNames',Animal);
    else
        AnimalTable = table(NumHypercolumns,MeanSpacing_mm,NumberPw,MeanPwDensity,RatioPrefCI,RatioSelectivityCI,RatioPwProb,RatioPwCI,RatioColumnSpacingCI,'RowNames',Animal);
    end
%     T = rows2vars(T);
%     T.Properties.VariableNames(1) = "Name";
%     AnimalTablePrint = rows2vars(AnimalTable);
%     AnimalTablePrint.Properties.VariableNames(1) = "Name";
% 
%     %% save table to latex
%     disp(AnimalTable)
%     table2latex(AnimalTable, [FigureFile 'Table.tex'])
    

end

function VectorArray = convert2vector(ArrayCell,data_obj,first)
    for ii = first:length(ArrayCell)
        vector = data_obj.array2vector(ArrayCell{ii});
        if ii == first
            VectorArray = zeros([length(vector) length(ArrayCell)]);
        end
        VectorArray(:,ii) = vector;
    end
end

function str = printValueCI(CIdown,Mean,CIup,precision)
    CImax = max(abs([CIdown-Mean CIup-Mean]));
    
    str = ['$' num2str(round(Mean,precision)) ' \pm ' num2str(round(CImax,precision)) '$'];
end


