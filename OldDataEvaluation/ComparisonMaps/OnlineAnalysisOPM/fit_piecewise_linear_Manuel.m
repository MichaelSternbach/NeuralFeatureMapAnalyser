function [pwd, plt, params, ssel, min_error] = fit_piecewise_linear_Manuel(x,y,fit_range,dplt,do_plotting)
% FUNCTION fit_piecewise_linear_WK_rerevised(x,y,do_plotting)
%  
%
% USAGE: [pwd, plt, params, ssel, min_error] = fit_piecewise_linear_WK_rerevised(x,y,do_plotting)
%
% Fits the stepwise linear function to y(x) so that it provides a best fit
% x, y are COLUMN vectors specifying the points to be fitted.
% The two vectors must be the same length.
%
% INPUT PARAMETERS:
% x             ... cutoff values in units of the local wavelength Lambda
% y             ... pinwheel density as a function of the cutoff wavelengths
% do_plotting   ... flag, if 1 plateau and fits are visualized
%
% OUTPUT PARAMETERS:
% pwd           ... pwd density, i.e. y-value of the fit at the plateau
% plt           ... plt ... [2x1] vector with beginning and end of the
%                   plateau
% params        ... actual fitting parameters
% ssel          ... fitting flag, 
%                   if [1 1 1], then there is a slope in the best fit on both
%                   sides of the plateau
%                   if [1 0 1], then there is a slope on right side of the
%                   plateau
%                   if [1 1 0], then there is a slope on left side of the
%                   plateau
% min_error     ... minimal squared error obtained by the fitting
%
%
%
%  Copyright (C) 2014, Max-Planck-Institut for Dynamics and Self-organization, The  
%  Nonlinear Dynamics Group. This software may be used, copied, or 
%  redistributed as long as it is not sold and this copyright notice is 
%  reproduced on each copy made. This routine is provided as is without 
%  any express or implied warranties whatsoever.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    min_error = 99999;
    
    if length(x) < 4
        disp('WTF');
    end
    
    mnx = fit_range(1); % min fit range
    mxx = fit_range(2); % max fit range
    %dplt = 0.4; %% minimal plateau length
    %%% These values are taken from Kaschube et al., Science 2010
    
    nexpn = 9;
    minplt = mnx;
    dxx = 0.05; % plateau increment
    dx = 0.05; % plateau shift increment
    nshif = 11; % maximal shift number
    
    errors = zeros(nshif, nexpn); 
    ppds = zeros(nshif, nexpn); 
    
    for ll  = 0:(nexpn-1) 
        dplat = dplt + ll *dxx;
        jj = 0;
        mxplt = minplt + dplat;
        
        while (mxplt < mxx && jj < nshif)
            
            sel = [1,1,1];%% Flag for the fitting
            mnplt = minplt + jj*dx;
            mxplt = minplt + jj*dx + dplat;
            
            if mnplt <= x(3)
                sel(2) = 0; %% Don't fit slope for values below plateau
            end
            if mxplt > x(end-3)
                sel(3) = 0; %% Don't fit slope for values above plateau
            end
            
            if sum(sel)>1
                %%% To the fitting with SVD
                [coeffs, mse, yy] = fitpiecewiselinearSVD( x', y',mnplt, mxplt, sel );
                errors(jj+1,ll+1) = mse;                
                ppds(jj+1,ll+1) = yy(find(x>mnplt+0.04,1));               
            end
            
            %%%% Store if error is minimal
            if mse < min_error
                min_error = mse;                
                params = coeffs;
                ssel = sel;
                jjmin = jj;
                llmin = ll;
                mnplt_min = mnplt;
                mxplt_min = mxplt;
            end
            jj = jj + 1;
        end
    end
    
    
    
    pwd = ppds(jjmin+1,llmin+1);
    plt = [mnplt_min mxplt_min];
    
    % To check how the best fit looks like
    if do_plotting == 1 && ~isdeployed       
        
        % Compute the actual fit function
        if isequal(ssel, [1, 1, 1])
            yy = params(1) + (-x+mnplt_min).*((-x+mnplt_min)>0)*params(2) + (-x+mxplt_min).*((-x+mxplt_min)<0)*params(3);       
        elseif isequal(ssel, [1, 0, 1])
            yy = params(1) +  (-x+mxplt_min).*((-x+mxplt_min)<0)*params(2);
        elseif isequal(ssel, [1, 1, 0])
            yy = params(1) + (-x+mnplt_min).*((-x+mnplt_min)>0)*params(2);
        end

        figure(101)
        plot(x,y,'x')
        hold on;
        plot(x,yy,'r')
        hold off;   
    
        pause(0.1);

    end



end
