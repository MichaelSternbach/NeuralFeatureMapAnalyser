function plotCI(CI,DoFilter,data_obj,folder,PwInfo,contourValueAngle,contourValueAbs)
    
    if DoFilter
        map = data_obj.filter_map(data_obj.read_map());
    else
        map = data_obj.read_map();
    end
    
    ROI=CI.BCA.ROI;
    
    if nargin <5
        PwInfo = 0;
    end
    
    if nargin < 7
        f = figure();
        tiledlayout(2,3)
        

        maxMap = 180;
        minMap = 0;
        plotAbsTile(CI.SE.CI_angle, 'CI angle SE [°]',maxMap,minMap,ROI,PwInfo)
        plotAbsTile(CI.BCA.CI_angle, 'CI angle BCA [°]',maxMap,minMap,ROI,PwInfo)


        nexttile;
        histogram(CI.BCA.CI_angle(ROI))
        xlabel('CI angle BCA [°]')
        ylabel('# number of pixels')

        [preMax,OrderMax] = getOrder(max([CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
        [preMin,OrderMin] = getOrder(max(-[CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
        maxMap = ceil(preMax)*10^OrderMax;
        minMap = -ceil(preMin)*10^OrderMin;
        plotAbsTile(CI.BCA.CI_Abs, 'CI abs BCA',maxMap,minMap,ROI,PwInfo)
        plotAbsTile(CI.SE.CI_Abs, 'CI abs SE',maxMap,minMap,ROI,PwInfo) 
        
        nexttile
        histogram(CI.BCA.CI_Abs(ROI))
        xlabel('CI abs BCA [°]')
        ylabel('# number of pixels')
        
        contourValueAngle = input('Input contour value for CI angle!');
        contourValueAbs = input('Input contour value for CI abs!');
        close(f)
        
    
    end
    %% final plot
    f = figure();
    tiledlayout(2,3)

    maxMap = 180;
    minMap = 0;
    plotAbsTile(CI.SE.CI_angle, 'CI angle SE [°]',maxMap,minMap,ROI,PwInfo)
    plotAbsTile(CI.BCA.CI_angle, 'CI angle BCA [°]',maxMap,minMap,ROI,PwInfo)
    
    
    nexttile;
    plotContourCI(contourValueAngle,map,ROI,CI.BCA.CI_angle,'angle','°')

    [preMax,OrderMax] = getOrder(max([CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
    [preMin,OrderMin] = getOrder(max(-[CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMin)*10^OrderMin;
    plotAbsTile(CI.BCA.CI_Abs, 'CI abs BCA',maxMap,minMap,ROI,PwInfo)
    plotAbsTile(CI.SE.CI_Abs, 'CI abs SE',maxMap,minMap,ROI,PwInfo)
    
    nexttile;
    plotContourCI(contourValueAbs,map,ROI,CI.BCA.CI_Abs,'abs','')

    if DoFilter
            CIFile = [folder 'Filtered_CI' data_obj.info.ID];
        else
            CIFile = [folder 'CI' data_obj.info.ID];
    end

    print(f,'-depsc', [CIFile '.eps'])
    savefig(f,[CIFile '.fig'])
end





