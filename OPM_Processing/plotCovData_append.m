function plotCovData_append(Covariances,data_obj,DoFilter,DiffType,FigureFile,MainTitle,mm,width_scale_pix)
    
    if nargin < 6
        MainTitle = {'Noise Covarinace Data', ' '};
        if DoFilter
            MainTitle = ['Filtered ' MainTitle];
        end
        if string(DiffType) == "align"
            MainTitle = ['Aligned ' MainTitle];
        end
    end
    if nargin <7
        mm = .5;
    end
    if nargin <8
        width_scale_pix = 5;
    end
    
    f = figure();
    t = tiledlayout(3,4);
    title(t,MainTitle)
    f.Position = [1 1 600 400];
    
    ROI = Covariances.ROI;
    spacing_pix = mm * data_obj.info.pix_per_mm*Covariances.scale;
    
    [YROI,XROI] = find(ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);
    
    %% plot variances
    [preMax,OrderMax] = getOrder(max([max(abs(Covariances.Var.C1),[],'all') max(abs(real(Covariances.Var.C2)),[],'all') max(abs(imag(Covariances.Var.C2)),[],'all')]));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMax)*10^OrderMax;%-ceil(preMin)*10^OrderMin;
    
    ax1 = nexttile;
    %[preMax,OrderMax] = getOrder(max(Covariances.Var.C1,[],'all'));
    plot_mapAbs(Covariances.Var.C1, 'Var C1',ceil(preMax)*10^OrderMax,-ceil(preMax)*10^OrderMax,ROI,ax1)
    hold on
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-red')
    hold on
    text(Xmin+width_scale_pix+spacing_pix,Ymin+width_scale_pix,[num2str(mm) ' mm'],'Color','red')
    xlim([Xmin, Xmax])
    ylim([Ymin, Ymax])
    
    
    ax1 = nexttile;
    %[preMax,OrderMax] = getOrder(max(Covariances.Var.C1,[],'all'));
    plot_mapAbs(imag(Covariances.Var.C1), 'Var Im C1',ceil(preMax)*10^OrderMax,-ceil(preMax)*10^OrderMax,ROI,ax1)
    hold on
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-red')
    hold on
    text(Xmin+width_scale_pix+spacing_pix,Ymin+width_scale_pix,[num2str(mm) ' mm'],'Color','red')
    xlim([Xmin, Xmax])
    ylim([Ymin, Ymax])
    
    %[preMax,OrderMax] = getOrder(max([real(Covariances.Var.C2) imag(Covariances.Var.C2)],[],'all'));
    %[preMin,OrderMin] = getOrder(max(-[real(Covariances.Var.C2) imag(Covariances.Var.C2)],[],'all'));
%     maxMap = ceil(preMax)*10^OrderMax;
%     minMap = -ceil(preMax)*10^OrderMax;%-ceil(preMin)*10^OrderMin;
    
    plotAbsTile(real(Covariances.Var.C2), 'Var Re C2',maxMap,minMap,ROI)
    hold on
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-red')
    hold on
    text(Xmin+width_scale_pix+spacing_pix,Ymin+width_scale_pix,[num2str(mm) ' mm'],'Color','red')
    xlim([Xmin, Xmax])
    ylim([Ymin, Ymax])
    
    plotAbsTile(imag(Covariances.Var.C2), 'Var imag C2',maxMap,minMap,ROI)
    hold on
    plot([Xmin+width_scale_pix,Xmin+width_scale_pix+spacing_pix],[Ymin+width_scale_pix, Ymin+width_scale_pix],'-red')
    hold on
    text(Xmin+width_scale_pix+spacing_pix,Ymin+width_scale_pix,[num2str(mm) ' mm'],'Color','red')
    xlim([Xmin, Xmax])
    ylim([Ymin, Ymax])

    %% plot 2d covariances
    [preMax,OrderMax] = getOrder(max([max(abs(Covariances.CoVar2D.C1),[],'all') max(abs(real(Covariances.CoVar2D.C2)),[],'all') max(abs(imag(Covariances.CoVar2D.C2)),[],'all')]));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMax)*10^OrderMax;%-ceil(preMin)*10^OrderMin;
    ROI = (Covariances.CoVar2D.N_PixelPairs>0);
    
%     [preMax,OrderMax] = getOrder(max(abs(Covariances.CoVar2D.C1),[],'all'));
    plotAbsTile(real(Covariances.CoVar2D.C1), 'Cov Re C1',maxMap,minMap,ROI)
    addScalBar(width_scale_pix,spacing_pix,mm)
    
    plotAbsTile(imag(Covariances.CoVar2D.C1), 'Cov Im C1',maxMap,minMap,ROI)
    addScalBar(width_scale_pix,spacing_pix,mm)
    
%     [preMax,OrderMax] = getOrder(max([real(Covariances.CoVar2D.C2) imag(Covariances.CoVar2D.C2)],[],'all'));
%     [preMin,OrderMin] = getOrder(max(-[real(Covariances.CoVar2D.C2) imag(Covariances.CoVar2D.C2)],[],'all'));
%     maxMap = ceil(preMax)*10^OrderMax;
%     minMap = -ceil(preMax)*10^OrderMax;%-ceil(preMin)*10^OrderMin;
    
    plotAbsTile(real(Covariances.CoVar2D.C2), 'Cov Re C2',maxMap,minMap,ROI)
    addScalBar(width_scale_pix,spacing_pix,mm)
    plotAbsTile(imag(Covariances.CoVar2D.C2), 'Cov Im C2',maxMap,minMap,ROI)
    addScalBar(width_scale_pix,spacing_pix,mm)
    
    ax2 = nexttile;
    MaxPixel = max(Covariances.CoVar2D.N_PixelPairs,[],'all');
    plot_mapAbs(Covariances.CoVar2D.N_PixelPairs, '# pixels',MaxPixel,0,ROI,ax2)
    addScalBar(width_scale_pix,spacing_pix,mm)
    
    %% plot 1D covariances
    
    NPixNonZero = find(movmean(Covariances.CoVar1D.NumDataList,5)==0);
    xMax = Covariances.CoVar1D.Distances(NPixNonZero(1))/data_obj.info.pix_per_mm/Covariances.scale;

    nexttile;
    plot(Covariances.CoVar1D.Distances/data_obj.info.pix_per_mm/Covariances.scale,Covariances.CoVar1D.C1,'DisplayName','C1')
    hold on
    plot(Covariances.CoVar1D.Distances/data_obj.info.pix_per_mm/Covariances.scale,real(Covariances.CoVar1D.C2),'DisplayName','Re C2')
    hold on
    plot(Covariances.CoVar1D.Distances/data_obj.info.pix_per_mm/Covariances.scale,imag(Covariances.CoVar1D.C2),'DisplayName','Im C2')
    xlabel('distance [mm]')
    xlim([0 xMax])
    legend()
    title('1D cov')

    nexttile;
    plot(Covariances.CoVar1D.Distances/data_obj.info.pix_per_mm/Covariances.scale,Covariances.CoVar1D.NumDataList)
    xlabel('distance [mm]')
    ylabel('# pixels')
    xlim([0 xMax])
    title('data for 1D cov')

%     if DoFilter
%             CovariancesFile = [folder DiffType 'Filtered_CoVar' data_obj.info.ID];
%         else
%             CovariancesFile = [folder DiffType 'CoVar' data_obj.info.ID];
%     end
    orient(f,'landscape')
    print(f, '-dpsc','-fillpage', '-append', [FigureFile '.ps'])
%     print(f,'-depsc', [CovariancesFile '.eps'])
%     savefig(f,[CovariancesFile '.fig'])
%     close

end

function addScalBar(width_scale_pix,spacing_pix,mm)
    hold on
    plot([width_scale_pix,width_scale_pix+spacing_pix],[width_scale_pix, width_scale_pix],'-red')
    hold on
    text(width_scale_pix+spacing_pix/2-8,width_scale_pix+4,[num2str(mm) ' mm'],'Color','red')
end