function plotPinwheelMovement(z,data_info,ROI,pinwheel_stats1,pinwheel_stats2,TrackingResults,FigureFileName)

    [YROI,XROI] = find(ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);

    f = figure();
    
    plot_map(z,ROI,0,1)
    
    hold on %pinwheel positions
    plot(pinwheel_stats1.x(:,1),pinwheel_stats1.y(:,1),'xw')
    % hold on;
    % contour(real(z),[0 0],'white')
    % hold on;
    % contour(imag(z),[0 0],'white')
    SizesCI = getConfidenceRegionPw(pinwheel_stats1,data_info.field_size_pix,0.95);
    
    %% plot pw IDs
    for ii = 1:size(TrackingResults.ini,1)
        if TrackingResults.ini(ii,2) ~= 0
            %text(pinwheel_stats.x(PwIDs(ii),1)+2,pinwheel_stats.y(PwIDs(ii),1),num2str(ii),'Color','white','FontSize',12)
            ID1 = TrackingResults.ini(ii,1);
            ID2 = TrackingResults.ini(ii,2);
            x0 = pinwheel_stats1.x(ID1,1);
            dx = pinwheel_stats2.x(ID2,1)-x0;
            y0 = pinwheel_stats1.y(ID1,1);
            dy = pinwheel_stats2.y(ID2,1)-y0;
            hold on
            quiver(x0, y0, dx, dy, 0, 'MaxHeadSize', 2, 'LineWidth', 2, 'Color', 'w');
        end
    end
    
    % % pinwheel density
    % title(['pinwheel density=' num2str(round(PwInfo.WeightedPwDensityFixedFilter,2)) '/\Lambda^2'])
    
    % save figure
    border = 5;
    xlim([Xmin-border Xmax+border])
    ylim([Ymin-border Ymax+border])

    print(f,'-depsc', FigureFileName)
end