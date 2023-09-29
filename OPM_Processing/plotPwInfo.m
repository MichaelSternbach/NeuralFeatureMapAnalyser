function plotPwInfo(PwInfo,data_obj,average_w,folder,MarkerSize)
    if nargin <5
        MarkerSize = 8;
    end
    
    if isfield(PwInfo,'CI_local_PwDensities')
        n_rows = 2;
        n_columns = 2;
    else
        n_rows = 1;
        n_columns = 3;
    end
    
    f = figure();
    tiledlayout(n_rows,n_columns)
    
    nexttile;
    plot_map(data_obj.filter_map(data_obj.read_map()),data_obj.ROI,0,1);
    hold
    plotPinwheelStats(PwInfo.pinwheel_stats,data_obj.info.field_size_pix)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    nexttile;
    [PWposY_plus, PWposX_plus] = find(PwInfo.pw_pos_matrix_plus > 0);
    [PWposY_minus, PWposX_minus]  = find(PwInfo.pw_pos_matrix_minus > 0);
    plot_map(data_obj.filter_map(data_obj.read_map()),data_obj.ROI,0,1);
%     set(gca, 'ydir','normal');
    hold on;
    plot(PWposX_plus,PWposY_plus,'o','color','w','MarkerSize',MarkerSize, 'linewidth', 2);
    plot(PWposX_minus, PWposY_minus,'o','color','k','MarkerSize',MarkerSize, 'linewidth', 2);
    axis image;
    hold off;

    nexttile;
    imagesc(PwInfo.automated_pw_density);    
%     set(gca, 'ydir','normal');
    axis image;
    colormap(jet);
    colorbar;
    hold on;
    plot(PWposX_plus,PWposY_plus,'o','color','w','MarkerSize',MarkerSize);
    plot(PWposX_minus, PWposY_minus,'o','color','k','MarkerSize',MarkerSize);
    hold off;
%     title('Automated pw dens estimates and pw positions');
    title('local pinwheel density');
    
    
    if isfield(PwInfo,'CI_local_PwDensities')
        
        nexttile;
        imagesc(abs(PwInfo.CI_local_PwDensities(:,:,1)-PwInfo.CI_local_PwDensities(:,:,2)));    
    %     set(gca, 'ydir','normal');
        axis image;
        colormap(jet);
        colorbar;
        hold on;
        plot(PWposX_plus,PWposY_plus,'o','color','w','MarkerSize',MarkerSize);
        plot(PWposX_minus, PWposY_minus,'o','color','k','MarkerSize',MarkerSize);
        hold off;
    %     title('Automated pw dens estimates and pw positions');
        title( ['CI pinwheel density for ' num2str(PwInfo.alpha) 'confidence']);
        
    end
%     nexttile;
%     histogram(PwInfo.d)
%     
%     nexttile;
%     histogram(PwInfo.d_eq)
%     
%     nexttile;
%     histogram(PwInfo.d_op)
    
    PwFigure = [folder 'PWs_' data_obj.info.ID];
    print(f,'-depsc', [PwFigure '.eps'])
    savefig(f,[PwFigure '.fig'])
        
    disp(['Pinwheel density with pw pos estimates is: ' num2str(length(PwInfo.PWxList)/(sum(data_obj.ROI(:))/(average_w*data_obj.info.pix_per_mm)^2))]);%understand!
    disp(['Pinwheel density with plateu fitting  is: ' num2str(PwInfo.pw_dens)])
        
    if isfield(PwInfo,'CI_PwDensities')
        disp(['The confidence interval for ' num2str(PwInfo.alpha) ' confidence is: ' num2str(abs(PwInfo.CI_PwDensities(1)-PwInfo.CI_PwDensities(2)))])
    end

end