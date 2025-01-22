function [x_range,y_range] = getRangeXY_ROI(ROI)

    [YROI,XROI] = find(ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);

    x_range = [Xmin Xmax];
    y_range = [Ymin Ymax];
end