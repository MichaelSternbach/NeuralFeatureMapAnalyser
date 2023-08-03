function [coeffs mse yy] = fitpiecewiselinearSVD( x, y,mnplt, mxplt, sel )
%
% USAGE: [coeffs mse yy] = fitpiecewiselinearSVD( x, y,mnplt, mxplt, sel )
%
% Fits the stepwise linear function to y(x) so that it provides a best fit
% produce a least-squares best fit to the data even if the data is
% overspecified or underspecified.
% much faster than the fminsearch procedure
% x, y are COLUMN vectors specifying the points to be fitted.
% The two vectors must be the same length.
%
% INPUT PARAMETERS:
% x        ... cutoff values in units of the local wavelength Lambda
% y        ... pinwheel density as a function of the cutoff wavelengths
% mnplt    ... minimum of plateau fitted (usually 0.2)
% mxplt    ... minimum of plateau fitted (usually 1.0)
% sel      ... flag for the fitting
%
% OUTPUT PARAMETERS:
% coeffs        ...
% ms            ...
% yy            ...
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    [sizexR, sizexC] = size(x);
    [sizeyR, sizeyC] = size(y);

    if (sizexC ~= 1)  || (sizeyC ~= 1)
        fprintf( 'Inputs of fit1dPolySVD must be column vectors\n' );
        return;
    end

    if (sizeyR ~= sizexR)
        fprintf( 'Inputs vectors of fit1dPolySVD must be the same length\n' );
        return;
    end


    numVals = sizexR;

    % Form array to process with SVD
    A = zeros(numVals,3);
    A(:,1) = x.*0+1;
    A(:,2) = (-x+mnplt).*((-x+mnplt)>0);
    A(:,3) = (-x+mxplt).*((-x+mxplt)<0);
    
    
    %%% Cuts the array, according to whether there is an upper or a lower
    %%% slope
       
    A = A(:,sel==1);
 

    % Perform SVD
    [u, s, v] = svd(A);

    % pseudo-inverse of diagonal matrix s
    sigma = eps; % minimum value considered non-zero
    qqs = diag(s);
    qqs(abs(qqs)>=sigma) = 1./qqs(abs(qqs)>=sigma);
    qqs(abs(qqs)<sigma)=0;
    qqs = diag(qqs);
    if numVals > 3
        qqs(numVals,1)=0; % add empty rows
    end

    % calculate solution
    coeffs = v*qqs'*u'*y;     
    
    % Calculate the error
    if isequal(sel, [1, 1, 1])
        yy = coeffs(1) + (-x+mnplt).*((-x+mnplt)>0)*coeffs(2) + (-x+mxplt).*((-x+mxplt)<0)*coeffs(3);       
    elseif isequal(sel, [1, 0, 1])
        yy = coeffs(1) +  (-x+mxplt).*((-x+mxplt)<0)*coeffs(2);
    elseif isequal(sel, [1, 1, 0])
        yy = coeffs(1) + (-x+mnplt).*((-x+mnplt)>0)*coeffs(2);
    else
        error('This function should be started with that sel flag!! ');
    end
    
    % Compute mean square error
    mse = sum((yy-y).^2)/length(yy);


end %%% END OF FUNCTION fitpiecewiselinearSVD.M

