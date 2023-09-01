function RectangleROI = getRectangleROI(rectangle,ROI)
    [X,Y]= meshgrid(1:size(ROI,2),1:size(ROI,1));
    RectangleROI = ROI.*(rectangle(3)>=X).*(X>=rectangle(1)).*(Y>=rectangle(2)).*(Y<=rectangle(4));
end