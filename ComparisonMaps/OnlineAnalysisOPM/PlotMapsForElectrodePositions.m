function PlotMapsForElectrodePositions(pinwheel_stats,orientation_stats,ElectrodePositions,data_obj,data_info,BloodVesselImg,data_path)

    %ElectrodePositions = {[border,border,border+lengthProngedElectrode_mm*pixels_per_mm,border]};
    %% Pinwheel Certainty positions
    FigureFilename = [data_path 'PW_Stats.fig'];
    figure
    z_base = data_obj.filter_map(data_obj.read_map());
    plot_map(z_base,data_obj.ROI,0,1)
    hold on
    plotPinwheelStats(pinwheel_stats,data_info.field_size_pix)
    PlotElectrodes(ElectrodePositions,'white')
    title('Pinwheel certainty')


    savefig([FigureFilename])

    %% Certainty OPM
    FigureFilename = [data_path 'OrientationStats'];
    figure
    tiledlayout(2,2)

    nexttile
    plot_map(orientation_stats(:,:,2),data_obj.ROI,0,1)
    hold on
    PlotElectrodes(ElectrodePositions,'white')
    title('orientation preference map')

    nexttile
    Abs = abs(orientation_stats(:,:,2))./mean(abs(orientation_stats(:,:,2)),'all');
    plot_mapAbs(Abs,'selectivity [<selectivity>]',max(Abs,[],'all'),min(Abs,[],'all'))
    PlotElectrodes(ElectrodePositions,'white')

    nexttile
    CI_Abs = abs(abs(orientation_stats(:,:,3))-abs(orientation_stats(:,:,1)))./mean(abs(orientation_stats(:,:,2)),'all');
    plot_mapAbs(CI_Abs,'uncertainty selectivity [<selectivity>]',max(CI_Abs,[],'all'),min(CI_Abs,[],'all'))
    PlotElectrodes(ElectrodePositions,'white')

    nexttile
    CI = abs(angle(orientation_stats(:,:,3)./orientation_stats(:,:,1)))/pi*90;
    plot_mapAbs(CI,'uncertainty orientation [°]',90,0)
    PlotElectrodes(ElectrodePositions,'white')
    
    savefig(FigureFilename)

    %% Certianty Contours
    figure
    tiledlayout(1,2)

    nexttile
    plotContourAngleDelta(20,orientation_stats(:,:,2),CI,data_obj.ROI)
    PlotElectrodes(ElectrodePositions,'white')
    
    nexttile
    plotContourAngleDelta(50,orientation_stats(:,:,2),CI,data_obj.ROI)
    PlotElectrodes(ElectrodePositions,'white')
    
    savefig([FigureFilename 'Contours'])

    %% Blood vessel maps
    figure
    tiledlayout(1,2)
    
    nexttile
    plot_map(z_base.*BloodVesselImg./abs(z_base.*BloodVesselImg),data_obj.ROI,0,1)
    hold on
    PlotBloodVessels(BloodVesselImg,data_obj.ROI,0.8)
    PlotElectrodes(ElectrodePositions,'red')
    
    nexttile
    PlotBloodVessels(BloodVesselImg,data_obj.ROI,1)
    PlotElectrodes(ElectrodePositions,'red')
    savefig([data_path 'BloodVessels+MapWithElectrodes'])
    
%     figure()
%     PlotBloodVessels(BloodVesselImg,ones(size(BloodVesselImg)))
%     PlotElectrodes(ElectrodePositions,'red')
%     savefig([data_path 'BloodVesselsWithElectrodes'])
end



function plotPinwheelStats(pinwheel_stats,field_size_pix)
    for i_pw = 1:getN_PW(pinwheel_stats)
        plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:),field_size_pix)
    end
end

function  plotPinwheel(PWx,PWy,ProbabilityPW,field_size_pix)
    ProbabilityLimitPW = .50;
    if ProbabilityPW >= ProbabilityLimitPW
        plotPosition(PWx(1),PWy(1),ProbabilityPW)
        plotConfidenceRegion(PWx,PWy,field_size_pix)
    end
end

function plotPosition(PWx,PWy,ProbabilityPW)
    plot(PWx,PWy,'.white')
    text(PWx,PWy,num2str(ProbabilityPW),'Color','white')
end
function plotConfidenceRegion(PWx,PWy,field_size_pix)
    CI = points_confidence_region(PWx,PWy,field_size_pix,'hull');
    contour(CI,[1 1],'white')
    %plotCI(CI)
end

function N_PW = getN_PW(pinwheel_stats)
    N_PW = size(pinwheel_stats.x,1);
end

function plotCI(CI)
    for i_x = 1:size(CI,1)
        for i_y = 1:size(CI,2)
            if CI(i_x,i_y)==1
                plot(i_x,i_y,'.white')
            end
        end
    end
end