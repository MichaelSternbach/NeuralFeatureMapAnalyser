function [PWNumbers,AreaSizes] = PwNuminCircl(z_input,ROI,RadiusLimits,NumPoints,seed,scale)
    rng(seed)
    
    disp('calculate pinwheels')
    [~,~,~,PWxList,PWyList,~, ~] = find_pinwheels(z_input,0,ROI);
    
    RandPoints = randi([1 sum(ROI,'all')],[NumPoints 1]);
    [ROIx,ROIy] = find(ROI);
    RandPointsX = ROIx(RandPoints);
    RandPointsY = ROIy(RandPoints);
    [X,Y]=meshgrid(1:size(ROI,2),1:size(ROI,1));
    
    disp('calculate pw numbers for random areas')
    PWNumbers = zeros(size(RandPoints));
    AreaSizes = zeros(size(RandPoints));
    %plot_map(ROI);hold on; plot(RandPointsX,RandPointsY,'+r')
    
%     p=10;
    for ii = 1:NumPoints
        Radius = (RadiusLimits(2)-RadiusLimits(1))*rand()+RadiusLimits(1);
        Distances = ((PWxList-RandPointsX(ii)).^2+(PWyList-RandPointsY(ii)).^2).^0.5;
        PWNumbers(ii) = sum(Distances<=Radius);
        %AreaSizes(ii) = sum(((ROIx-RandPointsX(ii)).^2+(ROIy-RandPointsY(ii)).^2).^0.5>=Radius,'all');%getAreaSize(RandPointsX(ii),RandPointsY(ii),ROI,Radius)
        Circle = (((X-RandPointsX(ii)).^2+(Y-RandPointsY(ii)).^2).^0.5>=Radius);
        AreaSizes(ii) = sum(ROI.*Circle,'all')*scale(RandPointsX(ii),RandPointsY(ii))^2;
        %plot_map(ROI+Circle.*1i);hold on; plot([RandPointsX(ii)],[RandPointsY(ii)],'+r')
        if AreaSizes(ii) == 0
            disp(['AreaSize = 0 for Radius = ' num2str(Radius)])
        end
%         if ii/size(NumPoints,2) >= p/100
%             disp([num2str(p),'%'])
%             p = p+10;
%         end
    end
end