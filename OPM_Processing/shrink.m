function A_shrink = shrink(A, shrinkage)
%
% FUNCTION CALLED BY FILTER_ORIENTATION_MAP.M
%
% USAGE : 
% A_shrink = shrink(A, shrinkage)
% This function shrinks a matrix by a certain factor shrink using local
% averaging of pixels
%
%
% INPUT ARGUMENTS:
% A .... two-dimensional matrix  to be shrunk
% shrinkage .... integer number larger than 0, 1 means no shrinkage at all,
% 2 means matrix is shrunk by a factor of 2, 
%
%
% OUTPUT ARGUMENTS: 
%
% A_shrink ... shrunk matrix
%
%  Copyright (C) 2014, Max-Planck-Institut for Dynamics and Self-organization, The  
%  Nonlinear Dynamics Group. This software may be used, copied, or 
%  redistributed as long as it is not sold and this copyright notice is 
%  reproduced on each copy made. This routine is provided as is without 
%  any express or implied warranties whatsoever.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    [n m] = size(A);

    % Cutoff pixles, such that mod(n, shrinkage == 0)
    while mod(n,shrinkage) ~= 0 
        A = A(1:n-1,:);
        n = n-1;

    end

    while mod(m,shrinkage) ~= 0 
        A = A(:,1:m-1);
        m = m-1;
    end

    A_shrink = zeros(n/shrinkage, m/shrinkage);
    
    new_i = 0;

    for i = 1 : shrinkage : n
        new_i = new_i + 1;
        new_j = 0;
        for j = 1 : shrinkage : m
            new_j = new_j + 1;        
            A_shrink(new_i, new_j) = sum(sum(A(i:i+shrinkage-1, j:j + shrinkage-1)))/shrinkage/shrinkage;        
        end
    end
    
    %%%%% Set the size of the map to an odd value (add -1's)
    A_shrink(:,size(A_shrink,2) + mod(size(A_shrink,2),2) ) = -1;
    A_shrink(size(A_shrink,1) + mod(size(A_shrink,1),2),: ) = -1;   

end



