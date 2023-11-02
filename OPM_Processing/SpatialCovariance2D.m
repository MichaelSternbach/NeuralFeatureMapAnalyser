function [C1,C2,N_PixelPairs] = SpatialCovariance2D(DiffMaps,ROI)

    InROI = reshape(ROI,[1 1 size(ROI)]).*reshape(ROI,[size(ROI) 1 1]);
    [InROIX,InROIY] = find(ROI);
    InROIX = reshape(InROIX,[1,length(InROIX)]);
    InROIY = reshape(InROIY,[1,length(InROIY)]);

    lengthROIx = max(InROIX,[],'all')-min(InROIX,[],'all');
    lengthROIy = max(InROIY,[],'all')-min(InROIY,[],'all');

    rangeX = -lengthROIx:lengthROIx;
    rangeY = -lengthROIy:lengthROIy;

    C1 = zeros(length(rangeX),length(rangeY));
    C2 = zeros(length(rangeX),length(rangeY));
    N_PixelPairs = zeros(length(rangeX),length(rangeY));




    disp('run 2D covariance calculation')
    p=10;
    for ix = 1:length(rangeX)
        for iy = 1:length(rangeY)
            [InROIX1,InROIY1,InROIX2,InROIY2]=getPixelPairsInROI(ROI,InROI,InROIX,InROIY,rangeX(ix),rangeY(iy));

            C1(ix,iy) = mean(DiffMaps(InROIX1,InROIY1,:).*conj(DiffMaps(InROIX2,InROIY2,:)),'all');
            C2(ix,iy) = mean(DiffMaps(InROIX1,InROIY1,:).*DiffMaps(InROIX2,InROIY2,:),'all');
            N_PixelPairs(ix,iy) = length(InROIX1);
        end
        %% disp progress
        if ix/length(rangeX) >= p/100
            disp([num2str(p),'%'])
            p = p+10;
        end
    end

%     figure
%     plot_mapAbs(real(C1))
% 
%     figure()
%     plot_mapAbs(real(C2))
% 
%     figure()
%     plot_mapAbs(imag(C2))


end

function [InROIX1,InROIY1,InROIX2,InROIY2]=getPixelPairsInROI(ROI,InROI,InROIX1,InROIY1,dX,dY)
    

    InROIX2 = InROIX1+dX;
    InROIY2 = InROIY1+dY;
    inFrame = find((0<InROIX2).*(InROIX2<=size(ROI,1)).*(0<InROIY2).*(InROIY2<=size(ROI,2)));
    inFrame = reshape(inFrame,[1 length(inFrame)]);

    InROIX1=InROIX1(inFrame);
    InROIY1=InROIY1(inFrame);
    InROIX2=InROIX2(inFrame);
    InROIY2=InROIY2(inFrame);
    
    if ~isempty(inFrame)

        SubPairs = sub2ind(size(InROI),InROIX1,InROIY1,InROIX2,InROIY2);
        PairsInROI = find(InROI(SubPairs));
        PairsInROI = reshape(PairsInROI,[1 length(PairsInROI)]);

        InROIX1=InROIX1(PairsInROI);
        InROIY1=InROIY1(PairsInROI);
        InROIX2=InROIX2(PairsInROI);
        InROIY2=InROIY2(PairsInROI);
    else
        InROIX1=[];
        InROIY1=[];
        InROIX2=[];
        InROIY2=[];
    end
    
%     %% plot to test
%     figure
%     plot_map(ROI)
%     hold on
%     plot(InROIX1,InROIY1,'+')
%     hold on
%     plot(InROIX2,InROIY2,'+')

%     InROI1 = sub2ind(size(ROI),InROIX1,InROIY1);
%     InROI2 = sub2ind(size(ROI),InROIX2,InROIY2);
end