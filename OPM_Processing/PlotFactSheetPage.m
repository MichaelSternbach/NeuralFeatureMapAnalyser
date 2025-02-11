function PlotFactSheetPage(animal,experiment_num,data_info,data_obj,DataFolder,FigureFile,...
    BiasDataFolder,width_scale_pix,color_contur,linewidth,Fontsize)


    
    %% set default values
    if nargin < 5
        FigureFile = [DataFolder 'Factsheet_' animal ' ' num2str(experiment_num) '.ps'];
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


    % if string(lower(animal)) == "dunnart"
    %     mm=0.2;
    % else
    %     mm=1;
    % end

    
    
    
    %% get column spacingmap and ID
    %[data_info,~,data_obj,~,BloodVesselImg] = getAnimalData(animal,experiment_num,AnimalDataFolder);
    [average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,DataFolder,false);
    z = data_obj.filter_map(data_obj.read_map());
    data_info.ID = replace(data_info.ID,'_',' ');


    %% load pinwheel data

    getCI = false;
    do_plotting=0;
    llp_cutoffs = linspace(0.01, 1,100);
    beta=0.5;

    Bootstrapsamples = 100;
    data_obj.prepare_samples_array(Bootstrapsamples)
    PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta);
    
    
    %% get CI data
    CI = calcCIs(data_obj,0.05,true,DataFolder);
    ROI = data_obj.ROI;


    %% make map borders ROI
    [YROI,XROI] = find(data_obj.ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);


    %% Borders Abs (CI) Plots
    [preMax,OrderMax] = getOrder(max([CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
    [preMin,OrderMin] = getOrder(max(-[CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMin)*10^OrderMin;


    %% make Figure
    close all
    f = figure();
    t = tiledlayout(5,3);
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
    
    

    plotAbsTile(CI.BCA.CI_angle, 'CI Prefernce [°]',maxMap,minMap,ROI)
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
    title('CI Prefernce CPDF')
    xlabel('CI Prefernce ≤ X [°]')
    ylabel('% of pix.')
    %xlim([0,1])
    axis(ax,'square')

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% plot polar map
    nexttile;
    scale = (data_info.pix_per_mm*average_spacing_mm).^-1;
    spacing_pix = 1/scale;
    plot_map(z,data_obj.ROI,1,1)
    hold on
    plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-white')
    hold on
    text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,'Λ','Color','white')
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color','white','linewidth',linewidth)
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
    title('Polar Map')

    %% plot selectivity Abs
    maxMap=7;
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
    title('rel. CI Selectivity CPDF')
    xlabel('CI Selectivity ≤ X [Selectivity]')
    ylabel('% of pix.')
    %xlim([0,1])
    axis(ax,'square')

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
    SizesCI = getConfidenceRegionPw(PwInfo.pinwheel_stats,data_info.field_size_pix,0.95);
    
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
    title('95% CI Pinwheel Positions')
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
    PwCI = SizesCI/(data_info.pix_per_mm)^2;
    plotCPDF(PwCI,'','-',ax)
    title('Pinwheel CI Size CPDF')
    xlabel('PW CI size ≤ X [mm^2]')
    ylabel('% of pinwheels')
    %xlim([0,1])
    axis(ax,'square')
    

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
    CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
    load(CIPwFile,'CI_PwDensities','alpha','PwInfosBS','PwInfosJS')
    LocalPwDensity = CI_PwDensities.LocalPwDensityPlateuFit(:,:,2);
    CI_LocalPwDensity = abs(CI_PwDensities.LocalPwDensityPlateuFit(:,:,3)-CI_PwDensities.LocalPwDensityPlateuFit(:,:,1));
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
    %[~,local_spacing_mm,~,~,~,~] =  loadColumnsSpacing(data_obj,DataFolder,true,true);

    ax = nexttile;
%         plot_mapAbs(local_spacing_mm,['Local Column Spacing [mm]'],max(local_spacing_mm(data_obj.ROI),[],'all'),min(local_spacing_mm(data_obj.ROI),[],'all'),data_obj.ROI,ax)
    plot_mapAbs(local_spacing_mm,['Local Column Spacing [mm]'],0.6,0.4,data_obj.ROI,ax)
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
    CISpacingFile = [DataFolder 'CI_MapSpacingFiltered_' data_obj.info.ID '.mat'];
    load(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm')
    
    AbsCI_local_spacing_mm = abs(CI_local_spacing_mm(:,:,2)-CI_local_spacing_mm(:,:,1));

    ax = nexttile;
    %plot_mapAbs(AbsCI_local_spacing_mm,['CI Local Column Spacing [mm]'],max(AbsCI_local_spacing_mm(data_obj.ROI),[],'all'),min(AbsCI_local_spacing_mm(data_obj.ROI),[],'all'),data_obj.ROI,ax)
    plot_mapAbs(AbsCI_local_spacing_mm,['CI Local Column Spacing [mm]'],0.5,0,data_obj.ROI,ax)
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
    title('rel. CI ColumnSpacing CPDF')
    xlabel('rel. CI \Lambda/ ≤ X ')
    ylabel('% of pix')
    %xlim([0,1])
    axis(ax,'square')
    
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
    RatioPwCI ={round(sum(PwCI<0.02,'all')/size(PwInfo.pinwheel_stats.probability,1)/size(PwInfo.pinwheel_stats.probability,2)*100,1)};
    RatioColumnSpacingCI = {round(sum(CI_CS<0.2,'all')/size(CI_CS,1)/size(CI_CS,2)*100,1)};
    
    if ~isempty(BiasDataFolder)
        load([BiasDataFolder 'PwDensityBiasCorrection.mat'],'PwDensityBiasCorrection')
        MeanPwDensityCorrected = {printValueCI(MeanPwDensityCI(1)-PwDensityBiasCorrection(experiment_num), PwInfo.WeightedPwDensityFixedFilter-PwDensityBiasCorrection(experiment_num), MeanPwDensityCI(2)-PwDensityBiasCorrection(experiment_num),2)};
        AnimalTable = table(MeanSpacing_mm,NumHypercolumns,NumberPw,MeanPwDensity,MeanPwDensityCorrected,RatioPrefCI,RatioSelectivityCI,RatioPwProb,RatioPwCI,RatioColumnSpacingCI,'RowNames',Animal);
    else
        AnimalTable = table(MeanSpacing_mm,NumberPw,MeanPwDensity,MeanPwDensityManuel,PlateauFitPwDensity,RatioPrefCI,RatioSelectivityCI,RatioPwProb,RatioPwCI,RatioColumnSpacingCI,'RowNames',Animal);
    end

    AnimalTable = rows2vars(AnimalTable);
    AnimalTable.Properties.VariableNames(1) = "Name";

    %% save table to latex
    disp(AnimalTable)
    table2latex(AnimalTable, [FigureFile 'Table.tex'])
    

    disp('Factsheet saved & finished')

end



function str = printValueCI(CIdown,Mean,CIup,precision)
    CImax = max(abs([CIdown-Mean CIup-Mean]));
    
    str = ['$' num2str(round(Mean,precision)) ' \pm ' num2str(round(CImax,precision)) '$'];
end



