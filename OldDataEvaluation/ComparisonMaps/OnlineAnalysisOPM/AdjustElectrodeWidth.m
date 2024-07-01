function pos = AdjustElectrodeWidth(positions,width)
    x1 = positions(1);
    y1 = positions(2);
    x2 = positions(3);
    y2 = positions(4);
    
    phi = atan((y2-y1)/(x2-x1));
    
    x2 = width*cos(phi)+x1;
    y2 = width*sin(phi)+y1;
    
    pos = [x1, y1, x2, y2];
end