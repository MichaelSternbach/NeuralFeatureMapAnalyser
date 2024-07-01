function PlotElectrodes(ElectrodePositions,color)
    for ii_Electrode = 1:size(ElectrodePositions,2)
        hold on
        ElectrodePosition = ElectrodePositions{ii_Electrode};
        if size(ElectrodePosition,2) == 2
            plot(ElectrodePosition(1),ElectrodePosition(2),'+','color',color)
            xEnd = ElectrodePosition(1);
            yEnd = ElectrodePosition(2);
        elseif size(ElectrodePosition,2) == 4
           plot(linspace(ElectrodePosition(1),ElectrodePosition(3),4),linspace(ElectrodePosition(2),ElectrodePosition(4),4),'-+','color',color) 
           xEnd = ElectrodePosition(3);
           yEnd = ElectrodePosition(4);
        end
        hold on
        if size(ElectrodePositions,2) > 1
            text(xEnd+1,yEnd+1,num2str(ii_Electrode),'Color',color)
        end
    end
end