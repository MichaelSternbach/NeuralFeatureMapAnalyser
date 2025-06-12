function TrackStats = plotPwTracks3D(PinwheelTracks,minSteps,ROI,plotID,alignTracks)

    if nargin <4
        plotID = false;
    end

    TrackStats.TrackSteps = zeros([1 size(PinwheelTracks.label,1)])*NaN;
    TrackStats.Start = zeros([1 size(PinwheelTracks.label,1)])*NaN;
    TrackStats.End = zeros([1 size(PinwheelTracks.label,1)])*NaN;
    TrackStats.TrackLength = zeros([1 size(PinwheelTracks.label,1)])*NaN;
    TrackStats.angle = zeros([1 size(PinwheelTracks.label,1)])*NaN;
    TrackStats.Tracks = cell([1 size(PinwheelTracks.label,1)]);
    disp(size(PinwheelTracks.label,1))
    n=0;
    for ii = 1:size(PinwheelTracks.label,1)
        x=PinwheelTracks.x(ii,:);
        y=PinwheelTracks.y(ii,:);
        z = 1:size(PinwheelTracks.label,2);
        
%         if nargin>2
%             insideROI = checkPathInROI(z, x, y, ROI);
%             % Set path points outside the ROI to NaN
%             x(~insideROI) = NaN;
%             y(~insideROI) = NaN;
%             z(~insideROI) = NaN;
%         end
        z(isnan(x)) = NaN;
        
        TrackStats.TrackSteps(ii) = sum(~isnan(x));
        idx=find(~isnan(x));
        if ~isempty(idx) && TrackStats.TrackSteps(ii)>minSteps
            TrackStats.Start(ii)=min(idx);
            TrackStats.End(ii)=max(idx);
            TrackStats.TrackLength(ii) = calculate3DPathLength(x,y,z);

            ang1 = angle((x(TrackStats.End(ii))-x(TrackStats.Start(ii)))+1i*(y(TrackStats.End(ii))-y(TrackStats.Start(ii))));
            TrackStats.angle1(ii) = ang1;

            r = sqrt((x(max(idx))-x(min(idx))).^2+(y(max(idx))-y(min(idx))).^2);
            dz = z(max(idx))-z(min(idx));
            ang2 = atan(r/dz);
            TrackStats.angle2(ii) = ang2;
        else
            TrackStats.Start(ii)=NaN;
            TrackStats.End(ii)=NaN;
            TrackStats.TrackLength(ii) = NaN;
            TrackStats.angle1(ii) = NaN;
            TrackStats.angle2(ii) = NaN;
        end

        if TrackStats.TrackSteps(ii)>minSteps
            if alignTracks
                plot3(x-x(min(idx)), y-y(min(idx)), z-z(min(idx)), 'LineWidth', 1);   % 'LineWidth' is optional for thicker lines
            else
                plot3(x, y, z, 'LineWidth', 1);   % 'LineWidth' is optional for thicker lines
            end
            n=n+1;
            hold on
        end
        
        x(isnan(x)) = [];
        y(isnan(y)) = [];
        z(isnan(z)) = [];

        if TrackStats.TrackSteps(ii)>minSteps && plotID
            text(x(1),y(1),z(1),[' ' num2str(ii)],"FontSize",12,'HorizontalAlignment', 'left')
            hold on
        end

        TrackStats.Tracks{ii}=permute([x;y;z],[2 1]);
           
    end
    disp(n)
    grid on; 

end

function insideROI = checkPathInROI(pathX, pathY, pathZ, roiMatrix)
    % CheckPathInROI determines whether points on a 3D path are inside a 3D ROI
    %
    % Inputs:
    %   pathX, pathY, pathZ - 1D vectors of the same length describing the 3D path
    %   roiMatrix - 3D boolean matrix representing the Region of Interest (ROI)
    %
    % Output:
    %   insideROI - boolean vector indicating if each point on the path is inside the ROI

    % Ensure the inputs are valid
    roiMatrix = (roiMatrix==1);
    if ~(isvector(pathX) && isvector(pathY) && isvector(pathZ))
        error('Path inputs must be 1D vectors.');
    end
    
    if length(pathX) ~= length(pathY) || length(pathY) ~= length(pathZ)
        error('Path vectors must have the same length.');
    end
    
    % Round the path coordinates to map them to ROI matrix indices
    pathX = round(pathX);
    pathY = round(pathY);
    pathZ = round(pathZ);

    pathX(pathX>size(roiMatrix,1)) = NaN;
    pathY(pathY>size(roiMatrix,2)) = NaN;
    pathZ(pathZ>size(roiMatrix,3)) = NaN;

% Get the dimensions of the ROI matrix
%     [roiSizeX, roiSizeY, roiSizeZ] = size(roiMatrix);
%     
%     % Check if each point is within the bounds of the ROI matrix
%     validIndices = (pathX >= 1 & pathX <= roiSizeX) & ...
%                    (pathY >= 1 & pathY <= roiSizeY) & ...
%                    (pathZ >= 1 & pathZ <= roiSizeZ);
%     
%     % Initialize the output vector
%     insideROI = false(size(pathX));
    
    % Only check valid indices within the bounds
%     for i = find(validIndices)'
    for i = 1:length(pathX)
        if isnan(pathX(i))||isnan(pathY(i))||isnan(pathZ(i))
            insideROI(i)=false;
        else
            insideROI(i) = roiMatrix(pathX(i), pathY(i), pathZ(i));
        end
    end
end


function l = calculate3DPathLength(x, y, z)
    % calculate3DPathLength Calculates the length of a 3D path
    % 
    % Inputs:
    %   x - Vector of x-coordinates
    %   y - Vector of y-coordinates
    %   z - Vector of z-coordinates
    %
    % Output:
    %   length - Total length of the 3D path

    % Validate inputs
    if length(x) ~= length(y) || length(y) ~= length(z)
        error('Input vectors x, y, and z must have the same length.');
    end

    % Calculate differences between consecutive points
    dx = diff(x);
    dy = diff(y);
    dz = diff(z);

    % Compute the Euclidean distances
    distances = sqrt(dx.^2 + dy.^2 + dz.^2);

    % Sum distances to get the total length
    l = nansum(distances);
end

