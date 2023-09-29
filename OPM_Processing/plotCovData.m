function plotCovData(Covariances,data_obj,DoFilter,DiffType,folder,MainTitle)
    
    if nargin < 6
        MainTitle = 'Covarinace data';
    end
    
    f = figure();
    t = tiledlayout(3,3);
    title(t,MainTitle)
    f.Position = [100 100 1100 1200];
    
    ROI = Covariances.ROI;
    ax1 = nexttile;
    
    [preMax,OrderMax] = getOrder(max(Covariances.Var.C1,[],'all'));
    plot_mapAbs(Covariances.Var.C1, 'Var C1',ceil(preMax)*10^OrderMax,0,ROI,ax1)
    
    
    [preMax,OrderMax] = getOrder(max([real(Covariances.Var.C2) imag(Covariances.Var.C2)],[],'all'));
    [preMin,OrderMin] = getOrder(max(-[real(Covariances.Var.C2) imag(Covariances.Var.C2)],[],'all'));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMin)*10^OrderMin;
    plotAbsTile(real(Covariances.Var.C2), 'Var Re C2',maxMap,minMap,ROI)
    plotAbsTile(imag(Covariances.Var.C2), 'Var imag C2',maxMap,minMap,ROI)


    
    ROI = (Covariances.CoVar2D.N_PixelPairs>0);
    
    [preMax,OrderMax] = getOrder(max(abs(Covariances.CoVar2D.C1),[],'all'));
    plotAbsTile(real(Covariances.CoVar2D.C1), 'Cov C1',ceil(preMax)*10^OrderMax,0,1)
    
    [preMax,OrderMax] = getOrder(max([real(Covariances.CoVar2D.C2) imag(Covariances.CoVar2D.C2)],[],'all'));
    [preMin,OrderMin] = getOrder(max(-[real(Covariances.CoVar2D.C2) imag(Covariances.CoVar2D.C2)],[],'all'));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMin)*10^OrderMin;
    plotAbsTile(real(Covariances.CoVar2D.C2), 'Cov Re C2',maxMap,minMap,ROI)
    plotAbsTile(imag(Covariances.CoVar2D.C2), 'Cov Im C2',maxMap,minMap,ROI)

    ax2 = nexttile;
    MaxPixel = max(Covariances.CoVar2D.N_PixelPairs,[],'all');
    plot_mapAbs(Covariances.CoVar2D.N_PixelPairs, '# pixels',MaxPixel,0,ROI,ax2)

    NPixNonZero = find(movmean(Covariances.CoVar1D.NumDataList,5)==0);
    xMax = Covariances.CoVar1D.Distances(NPixNonZero(1))/data_obj.info.pix_per_mm;

    nexttile;
    plot(Covariances.CoVar1D.Distances/data_obj.info.pix_per_mm,Covariances.CoVar1D.C1,'DisplayName','C1')
    hold on
    plot(Covariances.CoVar1D.Distances/data_obj.info.pix_per_mm,real(Covariances.CoVar1D.C2),'DisplayName','Re C2')
    hold on
    plot(Covariances.CoVar1D.Distances/data_obj.info.pix_per_mm,imag(Covariances.CoVar1D.C2),'DisplayName','Im C2')
    xlabel('distance [mm]')
    xlim([0 xMax])
    legend()
    title('1D cov')

    nexttile;
    plot(Covariances.CoVar1D.Distances/data_obj.info.pix_per_mm,Covariances.CoVar1D.NumDataList)
    xlabel('distance [mm]')
    ylabel('# pixels')
    xlim([0 xMax])
    title('data for 1D cov')

    if DoFilter
            CovariancesFile = [folder DiffType 'Filtered_CoVar' data_obj.info.ID];
        else
            CovariancesFile = [folder DiffType 'CoVar' data_obj.info.ID];
    end

    print(f,'-depsc', [CovariancesFile '.eps'])
    savefig(f,[CovariancesFile '.fig'])
%     close

end