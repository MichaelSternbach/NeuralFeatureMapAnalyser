function [NNDistances,NNDistances2S,NNDistancesOS] = getNNDistances(z_input,ROI,scale)
    if nargin ==2
        scale = ROI;
    end
    [~,~,~,PWxList,PWyList,signList, ~] = find_pinwheels(z_input,0,ROI);
    Signs = sign(signList);
    
    indxPW = sub2ind(size(scale),round(PWyList),round(PWxList));
    PwScales = scale(indxPW);
    
    Distances = ((PWxList-reshape(PWxList,[size(PWxList,2) 1])).^2+(PWyList-reshape(PWyList,[size(PWyList,2) 1])).^2).^0.5.*PwScales;
    SameSign = Signs.*reshape(Signs,[size(Signs,2) 1]) == 1;
    OppositeSign = Signs.*reshape(Signs,[size(Signs,2) 1]) == -1;


    Distances(find(Distances==0)) = inf;
    NNDistances = min(Distances,[],2);
    
    Distances2S = Distances;
    Distances2S(find((Distances.*SameSign)==0)) = inf;
    NNDistances2S = min(Distances2S,[],2);
    
    DistancesOS = Distances;
    DistancesOS(find((Distances.*OppositeSign)==0)) = inf;
    NNDistancesOS = min(DistancesOS,[],2);
end