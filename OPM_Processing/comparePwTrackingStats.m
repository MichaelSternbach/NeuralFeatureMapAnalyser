function PwTrackingStats = comparePwTrackingStats(pinwheel_stats1,pinwheel_stats2,SizesCI1,SizesCI2,TrackingResults)
    %% calc pw stats
    IDs1 = TrackingResults.ini(:,1);
    IDs2 = TrackingResults.ini(:,2);
    
    %% remove unmached pws
    IDs1 = IDs1(IDs2~=0);
    IDs2 = IDs2(IDs2~=0);
    
    %% calc distances and directions
    PwTrackingStats.distances = zeros(size(IDs1));
    PwTrackingStats.directions = zeros(size(IDs1));
    
    PwTrackingStats.prob1 = zeros(size(IDs1));
    PwTrackingStats.prob2 = zeros(size(IDs1));
    
    PwTrackingStats.SizeCI1 = zeros(size(IDs1));
    PwTrackingStats.SizeCI2 = zeros(size(IDs1));
    
    for ii = 1:length(IDs1)
        x1 = pinwheel_stats1.x(IDs1(ii),1);
        y1 = pinwheel_stats1.y(IDs1(ii),1);
    
        x2 = pinwheel_stats2.x(IDs2(ii),1);
        y2 = pinwheel_stats2.y(IDs2(ii),1);
    
        PwTrackingStats.distances(ii) = sqrt((x1-x2)^2+(y1-y2)^2);
        z = (x2-x1)+1i*(y2-y1);
        PwTrackingStats.directions(ii) = angle(z)/(2*pi)*360;
    
        PwTrackingStats.prob1(ii) = pinwheel_stats1.probability(IDs1(ii));
        PwTrackingStats.prob2(ii) = pinwheel_stats2.probability(IDs2(ii));
        
        PwTrackingStats.SizeCI1(ii) = SizesCI1(IDs1(ii));
        PwTrackingStats.SizeCI2(ii) = SizesCI2(IDs2(ii));
    end
end