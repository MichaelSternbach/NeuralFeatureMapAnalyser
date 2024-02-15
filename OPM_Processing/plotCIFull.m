function plotCIFull(CI,DoFilter,data_obj,PwInfo,folder,CovAligned,alpha)
       
    color_contur='w';
    linewidth = 0.1;
    
    %% set map & ROI
    z=data_obj.filter_map(data_obj.read_map());
    ROI=CI.BCA.ROI;
    [YROI,XROI] = find(ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);
    
    %% calc variances from aligned Cov complex field
    var_real = (real(CovAligned.Var.C1)-real(CovAligned.Var.C2))/2;
    var_imag = (real(CovAligned.Var.C1) + real(CovAligned.Var.C2))/2;
    
    Z = getZ(alpha);

    CI_Abs = (var_real).^0.5*Z*2;%/mean(abs(MeanMap),'all')

    CI_angle = asin((var_imag).^0.5)*Z*2;
    CI_angle(find(CI_angle>2*pi))=2*pi;
    CI_angle=real(CI_angle/pi*90);

    %% make function
    f = figure();
    s=1.5;
    f.Position = [1 1 800*s 200*s];
    tiledlayout(2,3)
    %% plot angle CI
    maxMap = 180;
    minMap = 0;
    plotAbsTile(CI.SE.CI_angle, 'CI angle SE [°]',maxMap,minMap,ROI,PwInfo,z,color_contur,linewidth,[Xmin Xmax],[Ymin Ymax])    
    plotAbsTile(CI.BCA.CI_angle, 'CI angle BCA [°]',maxMap,minMap,ROI,PwInfo,z,color_contur,linewidth,[Xmin Xmax],[Ymin Ymax])
    plotAbsTile(CI_angle, 'CI angle var [°]',maxMap,minMap,ROI,PwInfo,z,color_contur,linewidth,[Xmin Xmax],[Ymin Ymax])

    %% plot abs CI
    [preMax,OrderMax] = getOrder(max([CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
    [preMin,OrderMin] = getOrder(max(-[CI.BCA.CI_Abs CI.SE.CI_Abs],[],'all'));
    maxMap = ceil(preMax)*10^OrderMax;
    minMap = -ceil(preMin)*10^OrderMin;
    plotAbsTile(CI.BCA.CI_Abs, 'CI abs BCA',maxMap,minMap,ROI,PwInfo,z,color_contur,linewidth,[Xmin Xmax],[Ymin Ymax])
    plotAbsTile(CI.SE.CI_Abs, 'CI abs SE',maxMap,minMap,ROI,PwInfo,z,color_contur,linewidth,[Xmin Xmax],[Ymin Ymax])
    plotAbsTile(CI_Abs, 'CI abs var',maxMap,minMap,ROI,PwInfo,z,color_contur,linewidth,[Xmin Xmax],[Ymin Ymax])

    if DoFilter
            CIFile = [folder 'Filtered_CI_Full' data_obj.info.ID];
        else
            CIFile = [folder 'CI_Full' data_obj.info.ID];
    end

    print(f,'-depsc', [CIFile '.eps'])
    savefig(f,[CIFile '.fig'])
end





function Z = getZ(alpha) %https://www.mathsisfun.com/data/confidence-interval.html
    switch 1-alpha
        case .95
            Z = 1.960;
        case .99
            Z = 2.576;  
        otherwise
            disp(1-alpha)
    end
end