function [Distances,CoList,NumDataList] = SpatialCovariance(DiffMaps,ROI)
    InROI = reshape(ROI,[1 1 size(ROI)]).*reshape(ROI,[size(ROI) 1 1]);
    [X,Y] = meshgrid(1:size(DiffMaps,2),1:size(DiffMaps,1));
    X1 = reshape(X,[size(X) 1 1]);
    X2 = reshape(X,[1 1 size(X)]);
    Y1 = reshape(Y,[size(Y) 1 1]);
    Y2 = reshape(Y,[1 1 size(Y)]);
    Dsq = ((X1-X2).^2+(Y1-Y2).^2);
    SquaredDistances = reshape(reshape(0:size(ROI,1),[1 size(ROI,1)+1]).^2+reshape(0:size(ROI,2),[size(ROI,2)+1 1]).^2,[1 (size(ROI,1)+1)*(size(ROI,2)+1)]);
    SquaredDistances=unique(SquaredDistances.*(SquaredDistances<=max(Dsq,[],'all')));

    CoList = zeros(size(SquaredDistances));
    NumDataList = zeros(size(SquaredDistances));
    p=10;
    for ii = 1:size(SquaredDistances,2)
        arg = find((Dsq==SquaredDistances(ii)).*InROI);
        szDsq = size(Dsq);
        [I1,I2,I3,I4] = ind2sub(szDsq,arg);
        szDiff = size(DiffMaps);
        ArgDiffX = sub2ind(szDiff,repmat(I1,szDiff(3),1),repmat(I2,szDiff(3),1),repmat(reshape(1:szDiff(3),[szDiff(3) 1]),size(I2,1),1));
        ArgDiffY = sub2ind(szDiff,repmat(I3,szDiff(3),1),repmat(I4,szDiff(3),1),repmat(reshape(1:szDiff(3),[szDiff(3) 1]),size(I4,1),1));
        CoList(ii)=mean(DiffMaps(ArgDiffX).*conj(DiffMaps(ArgDiffY)),'all');
        NumDataList(ii) = size(ArgDiffX,1);
        %% disp progress
        if ii/size(SquaredDistances,2) >= p/100
            disp([num2str(p),'%'])
            p = p+10;
        end
    end
    Distances = SquaredDistances.^0.5;
end