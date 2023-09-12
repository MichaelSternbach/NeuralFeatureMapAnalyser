function  [all_areas, all_n, areas ] = gimme_nv_roi_with_local_pw_dens(pw_pos_matrix,ITER, lambda, roi_matrix)
% DESCRIPTION: Computes number variance, standard deviation and count variance as a function of circular area size
%               using the local pinwheel density matrix and not the
%               estimated pinwheel positions
%
%
% INPUT PARAMETERS:
%
% pw_pos_matrix ... local piwnheel density of map, computed with 
% sigma = 0.01Lambda used to for to determine the number of pinwheels in a given area 
% ITER ... number of circles used for each radius
% lambda in units of pixels
% roi_matrix ... is one inside the ROI and zero otherwise
%
% OUTPUT PARAMETERS
% all_areas ... overlaps of ROI and circle
% all_n ... ... number of pinwheels within that overlap, computed with
% local pinwheel density determined with sigma = 0.01
% 
%
%

   
    rand('seed', 1234);    
    
    [X, Y] = meshgrid(1:size(roi_matrix,2),1:size(roi_matrix,1));

    
    % Determine extension of the roi to draw the positions of the circles    
    [x, y] = find(roi_matrix == 1);
    x_roi_min = min(x);
    x_roi_max = max(x);
    y_roi_min = min(y);
    y_roi_max = max(y);
    
    
    %%%%%%%%% Draw positions of circles centers
    x_centers = rand(1,ITER)*(x_roi_max - x_roi_min) + x_roi_min;
	y_centers = rand(1,ITER)*(y_roi_max - y_roi_min) + y_roi_min;
	
    
    s = size(roi_matrix,1);
    
    %%% Radii in units of lambda
    %%% Maximum radius is half the size of the matrix
	ll = logspace(-1,log10(min(100,s/lambda/2)),200);
    %ll = sqrt(linspace(.1^2,min(100,s/l_typ/2)^2,500));
    
    all_areas = [];
    all_n = [];
    
    disp('Starting the estimate count variance for nonperiodic OPM...');
    
    for l = ll %% radius in intr. units
		
        % Area is returned in pixel size
        % Radii should be in units of pixels
        [n, areas] = how_many_dist_2d(pw_pos_matrix,l*lambda,x_centers,y_centers, X,Y, roi_matrix);
        
        % Add only those, for which there has been a real overlap and a pinwheel
        all_areas = [all_areas areas(areas>0)/lambda/lambda];
        all_n = [all_n n(areas>0)];
	            
    end
    disp('...done.')
    
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%   LIST OF USED FUNCTIONS            %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n areas]  = how_many_dist_2d(pw_pos_matrix,radius,x_centers,y_centers, X,Y, roi_matrix)

	ITER = length(x_centers);

	
    areas = zeros(1, ITER);
    n = zeros(1, ITER);
	
	for t = 1:ITER
        
        [areas(t) overlap_mat] = compute_overlap_circle_ROI(roi_matrix,X,Y, x_centers(t), y_centers(t),radius);       
        n(t) = sum(pw_pos_matrix(overlap_mat(:) == 2));

	end
end



function [overlap overlap_mat] = compute_overlap_circle_ROI(roi_matrix,X,Y, x_cent, y_cent,radius)
% radius and xp, yp should be given in pixel units

    
    overlap_mat = roi_matrix + (sqrt((X-x_cent).^2 + (Y-y_cent).^2) < radius);    
    overlap = sum(overlap_mat(:) == 2);
end
% 
% function circle(x,y,r)
%     %x and y are the coordinates of the center of the circle
%     %r is the radius of the circle
%     %0.01 is the angle step, bigger values will draw the circle faster but
%     %you might notice imperfections (not very smooth)
%     ang=0:0.01:2*pi; 
%     xp=r*cos(ang);
%     yp=r*sin(ang);
%     plot(x+xp,y+yp, 'g', 'Linewidth', 2);
% end