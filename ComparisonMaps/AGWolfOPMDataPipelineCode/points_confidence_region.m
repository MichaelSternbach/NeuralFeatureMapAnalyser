function CI = points_confidence_region(X,Y,size_area,method)
% methods: 'hull' or 'gaussians'

if numel(size_area)==1
    size_area = [size_area,size_area];
end

if nargin<4
    method = 'hull';
end

% clean
X=X(~isnan(X));
Y=Y(~isnan(Y));

CI = false(size_area);

switch method
    case 'gaussians'
        
        [Xr,Yr] = meshgrid(1:size_area(1),1:size_area(2));
        sig = 1;
        pw_dens = false(size_area);
        for ii = 1:length(X)
            pw_dens = pw_dens + exp(-1/2/sig^2*((Xr-X(ii)).^2+(Yr-Y(ii)).^2));
        end
        
        [val,ind]=sort(pw_dens(:)/sum(pw_dens(:)),'descend');
        th = find(cumsum(val)>0.95,1,'first');
        CI(ind(1:th)) = true;
        
        % fill holes and pick largest area
        CI = imfill(CI,'holes');
        CC = bwconncomp(CI);
        [~,ind]=max(cellfun(@length,CC.PixelIdxList));
        CI = false(size_area);
        CI(CC.PixelIdxList{ind})=true;
    
    case 'boot_gaussians'
        
        [Xr,Yr] = meshgrid(1:size_area(1),1:size_area(2));
        sig = 1;
        pw_dens = false(size_area);
        for ii = 1:1000
            ind = ceil(length(X)*rand);
            pw_dens = pw_dens + exp(-1/2/sig^2*((Xr-X(ind)).^2+(Yr-Y(ind)).^2));
        end
        
        [val,ind]=sort(pw_dens(:)/sum(pw_dens(:)),'descend');
        th = find(cumsum(val)>0.95,1,'first');
        CI(ind(1:th)) = true;
        
        % fill holes and pick largest area
        CI = imfill(CI,'holes');
        CC = bwconncomp(CI);
        [~,ind]=max(cellfun(@length,CC.PixelIdxList));
        CI = false(size_area);
        CI(CC.PixelIdxList{ind})=true;
        
    case 'hull'
        % make a hull and cover all points
        try
            if length(unique(X))>=3 || length(unique(X))>=3
                
                K = convhull(X,Y);
                % include center of hull
                CI = poly2mask(X(K),Y(K),size_area(1),size_area(2));
                % include hull lines
                for ki = 1:length(K)-1
                    [x,y]=bresenham(X(K(ki)),Y(K(ki)),X(K(ki+1)),Y(K(ki+1)));
                    CI(sub2ind(size(CI),y,x))=true;
                end
                [x,y]=bresenham(X(K(end)),Y(K(end)),X(K(1)),Y(K(1)));
                CI(sub2ind(size(CI),y,x))=true;
            end
        catch
            % disp('Not enough  points to make hull')
        end
        
        % include the pixel where the points belong to
        CI(sub2ind(size(CI),ceil(Y(:)),ceil(X(:)))) = true;
        CI(sub2ind(size(CI),floor(Y(:)),floor(X(:)))) = true;
        CI(sub2ind(size(CI),floor(Y(:)),ceil(X(:)))) = true;
        CI(sub2ind(size(CI),ceil(Y(:)),floor(X(:)))) = true;
        
end

end


function [x,y]=bresenham(x1,y1,x2,y2)
%Matlab optmized version of Bresenham line algorithm. No loops.
%Format:
%               [x y]=bham(x1,y1,x2,y2)
%
%Input:
%               (x1,y1): Start position
%               (x2,y2): End position
%
%Output:
%               x y: the line coordinates from (x1,y1) to (x2,y2)
%
%Usage example:
%               [x y]=bham(1,1, 10,-5);
%               plot(x,y,'or');
x1=round(x1); x2=round(x2);
y1=round(y1); y2=round(y2);
dx=abs(x2-x1);
dy=abs(y2-y1);
steep=abs(dy)>abs(dx);
if steep
    t=dx;
    dx=dy;
    dy=t;
end

%The main algorithm goes here.
if dy==0
    q=zeros(dx+1,1);
else
    q=[0;diff(mod([floor(dx/2):-dy:-dy*dx+floor(dx/2)]',dx))>=0];
end

%and ends here.

if steep
    if y1<=y2
        y=[y1:y2]';
    else
        y=[y1:-1:y2]';
    end
    
    if x1<=x2
        x=x1+cumsum(q);
    else
        x=x1-cumsum(q);
    end
else
    if x1<=x2
        x=[x1:x2]';
    else
        x=[x1:-1:x2]';
    end
    if y1<=y2
        y=y1+cumsum(q);
    else
        y=y1-cumsum(q);
    end
end
end

