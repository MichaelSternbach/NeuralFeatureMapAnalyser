function PlotPwCI(animal,data_info,data_obj,DataFolder)
    close all
    
    rectangle = false;
    
    
    if string(lower(animal)) == "dunnart"
        mm=0.2;
    else
        mm=1;
    end

    width_scale_pix= 15;

    
    %% figure file name
    FigureFile = [DataFolder 'Factsheet_' animal ' ' data_info.ID]
    
    
    %% contour plots
    grey = [0.6 0.6 0.6];
    color_contur = grey;%'grey';
    linewidth = .3;
    set(gca,'Fontsize',20)


    %% animal
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
    t = tiledlayout(5,3);
    s=1;
    f.Position = [100 100 594*s 841*s];
    title(t,data_info.ID)


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
    
    


    %% Borders Abs (CI) Plots
    [preMax,OrderMax] = getOrder(max([CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
    [preMin,OrderMin] = getOrder(max(-[CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMin)*10^OrderMin;
    
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


    xlim([Xmin Xmax])
    ylim([Ymin Ymax])

    yticks([])
    xticks([])
    
    %% pw prob. and confidence region

    ax = nexttile;
    imagesc(ax,zeros(size(ROI)))
    hold on;
    SizesCI = getConfidenceRegionPw(PwInfo.pinwheel_stats,data_info.field_size_pix,0.95);
    
    contour(ROI,[1 1],'white','linewidth',linewidth)
    m=100;
    cm_viridis=viridis(m);
    colormap(ax, cm_viridis);
    axis image
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
    title('95% CI Pinwheel Positions')
    yticks([])
    xticks([])
           
    
            
    %% plot CPDF pinwheel CI Size
    ax = nexttile;
    PwCI = SizesCI/(data_info.pix_per_mm)^2;
    plotCPDF(PwCI,'','-',ax)
    title('Pinwheel CI Size CPDF')
    xlabel('PW CI size ≤ X [mm^2]')
    ylabel('% of pinwheels')
    %xlim([0,1])
    axis(ax,'square')
    

%     %% load pinwheel CI data
%     CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
%     load(CIPwFile,'CI_PwDensities','alpha','PwInfosBS','PwInfosJS')
%     LocalPwDensity = CI_PwDensities.LocalPwDensityPlateuFit(:,:,2);
%     CI_LocalPwDensity = abs(CI_PwDensities.LocalPwDensityPlateuFit(:,:,3)-CI_PwDensities.LocalPwDensityPlateuFit(:,:,1));
%     
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
    
    

    
    %% plot CPDF pinwheel CI Size
    ax = nexttile;
    CI_CS = AbsCI_local_spacing_mm./(local_spacing_mm*2);
    plotCPDF(CI_CS(ROI),'','-',ax)
    title('rel. CI ColumnSpacing CPDF')
    xlabel('rel. CI \Lambda/ ≤ X ')
    ylabel('% of pix')
    %xlim([0,1])
    axis(ax,'square')
    
    print(f, '-depsc', [FigureFile '.eps'])      
    
end



function str = printValueCI(CIdown,Mean,CIup,precision)
    CImax = max(abs([CIdown-Mean CIup-Mean]));
    
    str = ['$' num2str(round(Mean,precision)) ' \pm ' num2str(round(CImax,precision)) '$'];
end



