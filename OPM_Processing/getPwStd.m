function PwStd = getPwStd(pinwheel_stats)
    PwStd = sqrt(nansum((pinwheel_stats.x(:,1)-pinwheel_stats.x).^2+(pinwheel_stats.y(:,1)-pinwheel_stats.y).^2,2)./(sum(~isnan(pinwheel_stats.x),2)-1));
    
end