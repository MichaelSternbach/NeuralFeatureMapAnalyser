function plotGrayMap(map)
    figure()
    a = abs(map);
    %figure;
    imagesc(a); 
    %colormap jet;
    %colormap hsv;
    colormap gray;
end