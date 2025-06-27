function PlotFiguresCI(animal, spec_num,AnimalDataFolder,DataFolder,FigureFolder,width_scale_pix,color_contur,linewidth,Fontsize)
    
    %% define data folder
    DataFolderAnimal = [DataFolder animal '/' animal num2str(spec_num) '/'];

    %% default FigureFolder
    if nargin <5
        FigureFolder = DataFolderAnimal;
    end

    %% defaul plot parameter
    if nargin < 8
        color_contur = [0.6 0.6 0.6];%%'grey';
    end
    if nargin < 9
        linewidth = 0.5;
    end
    % if nargin < 10
    %     Fontsize = 20;
    % end
    % if nargin < 7
    %     width_scale_pix = 15;
    % end

    %% get animal data
    [~,~,data_obj,~,~] = getAnimalData(animal,spec_num,AnimalDataFolder);
    ROI = data_obj.ROI;
    data_obj.activateGIF(true,4)
    z = data_obj.filter_map(data_obj.read_map());


    %% find map borders ROI
    [YROI,XROI] = find(data_obj.ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);

    %% get column spacing
    SpacingFile = [DataFolderAnimal 'MapSpacingFiltered_' data_obj.info.ID '.mat'];
    load(SpacingFile,'average_spacing_mm')


    %% load CI data   
    CI = calcCIs(data_obj,0.05,true,DataFolderAnimal);

    PwInfoFile = [DataFolderAnimal 'PwInfo_' data_obj.info.ID '.mat'];
    load(PwInfoFile,'PwInfo')

    % CIPwFile = [DataFolderAnimal 'CI_PwDensity_' data_obj.info.ID '.mat'];
    % load(CIPwFile,'CI_PwDensities','alpha','PwInfosBS','PwInfosJS')
    
    % CISpacingFile = [DataFolderAnimal 'CI_MapSpacingFiltered_' data_obj.info.ID '.mat'];
    % load(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm')

    
    %% plot CI preference
    figure;
    plotAbsTile(CI.BCA.CI_angle, 'CI Prefernce [°]',180,0,ROI)
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
    hold on
    contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    hold on
    contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)

    xlim([Xmin Xmax])
    ylim([Ymin Ymax])

    print('-depsc', [FigureFolder data_obj.info.ID 'PreferenceCI.eps'])
    
    


    %% plot CI Selectivity
    figure;
    CI_Abs = CI.BCA.CI_Abs;
    maxMap = max(CI_Abs(ROI),[],'all');
    minMap = min(CI_Abs(ROI),[],'all');
    plotAbsTile(CI_Abs, 'CI abs BCA',maxMap,minMap,ROI)
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
    hold on
    contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    hold on
    contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)

    xlim([Xmin Xmax])
    ylim([Ymin Ymax])

    print('-depsc', [FigureFolder data_obj.info.ID 'SelectivityCI.eps'])

    %% plot Pw CI

    figure;
    sigma = 0.05;
    PWxList = PwInfo.pinwheel_stats.x(~isnan(PwInfo.pinwheel_stats.x));
    PWyList = PwInfo.pinwheel_stats.y(~isnan(PwInfo.pinwheel_stats.y));

    local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PWxList, PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
    local_pw_dens = local_pw_dens./sum(local_pw_dens,'all').*sum(PwInfo.pinwheel_stats.probability,'all');

    plot_mapAbs(local_pw_dens,'Pinwheel Prob. Density',max(local_pw_dens(data_obj.ROI),[],'all'),min(local_pw_dens(data_obj.ROI),[],'all'),data_obj.ROI)
    
    hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'.','Color',color_contur,'linewidth',linewidth)
    hold on
    contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    hold on
    contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)

    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
    yticks([])
    xticks([])

    print('-depsc', [FigureFolder data_obj.info.ID 'PwDensityCI.eps'])
        
end