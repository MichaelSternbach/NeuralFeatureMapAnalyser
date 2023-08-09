function [PWDensitySD,NumAreas] = getPWDensitySD(PWNumbers,AreaSizes,AreaList,dA)
    PWDensitySD = zeros(size(AreaList));
    NumAreas = zeros(size(AreaList));
    for ii = 1:size(AreaList,2)
        ArgPWNum = find((AreaSizes>=AreaList(ii)).*(AreaSizes<=(AreaList(ii)+dA)));
        PWDensitySD(ii)=std(PWNumbers(ArgPWNum),0,'all');
        NumAreas(ii) = sum(ArgPWNum,'all');
    end
end