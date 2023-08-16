function [NNDistances,NNDistances2S,NNDistancesOS] = NNDistances(z_input,ROI)
    [~,~,~,PWxList,PWyList,signList, ~] = find_pinwheels(z_input,0,ROI);
    Signs = sign(signList);

    Distances = ((PWxList-reshape(PWxList,[size(PWxList,2) 1])).^2+(PWyList-reshape(PWyList,[size(PWyList,2) 1])).^2).^0.5;
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