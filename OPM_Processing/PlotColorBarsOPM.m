close all

n = 8;

figure
% for i = 0:n
%     makeRectangle(0,0,180/n*i)
%     hold on
% end

makeRectangle(1,4,0)
hold on
makeRectangle(1,0,0)
axis equal; axis off

print('-depsc', 'OPM_Colorbar.eps')
close






function makeRectangle(x,y,angle)
angle = 2*pi*angle/360;

lx = 4;
ly = 1;

r = sqrt((lx/2)^2+(ly/2)^2);
angleRectangle = acos(lx/2/r);

h = rectangle('Position',[x y x+lx y+ly],'FaceColor',getRGB_colorsOPM(2*angle));%[x-cos(angleRectangle)*r y-sin(angleRectangle)*r x+cos(angleRectangle)*r y+sin(angleRectangle)*r]
%g = hgtransform('Matrix', makehgtform('zrotate', angle));
%Sh.Parent = g;
end

function RGB = getRGB_colorsOPM(angle)

cm = makeColormap('orientation',16);

ori=(angle)/(2*pi);ori(ori<0)=1+ori(ori<0);

% find to which color interval it belongs
intervals=linspace(0,1,2*size(cm,1)+1);
[~,ori] = histc(ori,intervals(2:2:end));
ori = ori+1;
ori(ori==size(cm,1)+1)=1;

RGB = cm(ori(:),:);%[1 0 0];
end