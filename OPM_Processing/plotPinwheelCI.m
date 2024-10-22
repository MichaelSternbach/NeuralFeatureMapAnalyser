function SizesCI = plotPinwheelCI(z,data_info,ROI,pinwheel_stats,PwIDs,FigureFileName)

    [YROI,XROI] = find(ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);

    f = figure();
    
    plot_map(z,ROI,0,1)
    
    hold on %pinwheel positions
    plot(pinwheel_stats.x(:,1),pinwheel_stats.y(:,1),'xw')
    % hold on;
    % contour(real(z),[0 0],'white')
    % hold on;
    % contour(imag(z),[0 0],'white')
    
    %% plot position confidence intervall
    hold on;
    SizesCI = getConfidenceRegionPw(pinwheel_stats,data_info.field_size_pix,0.95);
    
    %% plot pw IDs
    for ii = 1:size(PwIDs)
        if PwIDs(ii) ~= 0
            text(pinwheel_stats.x(PwIDs(ii),1)+2,pinwheel_stats.y(PwIDs(ii),1),num2str(ii),'Color','white','FontSize',12)
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