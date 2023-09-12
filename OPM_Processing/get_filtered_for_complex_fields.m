function [map_dataf,filt] = get_filtered_for_complex_fields(map,roi_matrix,filter_type,cut_off,beta)
%
% FUNCTION CALLED BY FILTER_ORIENTATION_MAP.M AND ESTIMATE_LOCAL_PW_DENSITY.M
%
% USAGE : 
% 
% GET_FILTERED_FOR_COMPLEX_FIELDS(map_data,roi_matrix,llp,hhp,beta)
%
%
% GET_FILTERED_FOR_COMPLEX_FIELDS(map_data,roi_matrix,llp,hhp)
% works with GAUSSIAN filters
% Returns filtered image and filter kernel. llp is the lower frequency
% cutoff, hhp the high-frequency cutoff
%
% GET_FILTERED_FOR_COMPLEX_FIELDS(map_data,roi_matrix,llp,hhp, beta)
% will work with with FERMI filters
% Returns filtered image and filter kernel. 
%
%
% INPUT PARAMETERS:
% map ...      complex-valued two-dimensional array
% roi_matrix ...    matrix consisting of 1 within the region of interest, 0
%                   otherwise
% filter_type   ... string, choose either 'lowpass' or 'highpass'
% cut_off       ... cut off of filter in millimeter
% beta          ... beta specifies the slope of the fermi function if
%                   Fermi-Filter is used
% box_filter    ... flag, if set to 1, Fermi-Filter is thresholded at 0.5
%                   to obtain a box filter
%
% OUTPUT PARAMETERS: 
%   map_dataf   ... filtered map
%   filt        ... filter used for filtering (in Fourier Space)
%
%
%
%  Copyright (C) 2014, Max-Planck-Institut for Dynamics and Self-organization, The  
%  Nonlinear Dynamics Group. This software may be used, copied, or 
%  redistributed as long as it is not sold and this copyright notice is 
%  reproduced on each copy made. This routine is provided as is without 
%  any express or implied warranties whatsoever.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ext = round(2.^(ceil(log(max(size(map)))/log(2)))); %(This will be the size of the filters/arrays)
   
    a = size(map);

    bg=zeros(ext);%put on square array of size 2^N (zero padding)
    bg(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)))=map;
    roi_padd = zeros(ext);
    roi_padd(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)))=roi_matrix;

    roi = find(roi_padd == 1);
    [xx,yy] = meshgrid(1:ext,1:ext);%Koordinate system
    xx = xx-ext/2-1;
    yy = yy-ext/2-1;
 
    filt = zeros(size(bg)); 
    
    
    if cut_off > 0
    
        if nargin > 4 %% --> implies Fermi Filtering

            cut_off=ext./cut_off;
            beta = beta*cut_off;
            filt = 1./(1+exp((sqrt(xx.^2+yy.^2)-cut_off)./(beta)));
            
            if nargin > 5 %% ---> implies that box_off options is specified
                if box_filter == 1
                    disp('Using Fermi-derived box filter!');
                    indq = filt > 0.5;
                    filt = filt*0;
                    filt(indq) = 1;
                end
            end
            
        else %% --> implies Gaussian Filtering
            cut_off = ext/cut_off;
            filt = exp(-(xx.^2+yy.^2)/2/cut_off^2);
        end


        map_dataf=fftshift(ifft2(fft2(fftshift(bg)).*fftshift(filt)));
        %%% Take Borders into consideration
        overlap = fftshift(ifft2(fft2(fftshift(roi_padd)).*fftshift(filt)));       
        if strcmpi(filter_type, 'lowpass') || strcmpi(filter_type, 'low_pass') % Fermi Lowpass
            map_dataf(roi) = map_dataf(roi)./overlap(roi);
        elseif strcmpi(filter_type, 'highpass') || strcmpi(filter_type, 'high_pass') % Fermi Highpass
            map_dataf(roi) = bg(roi) -  map_dataf(roi)./overlap(roi);
        else
            error('Invalid filter type: Choose either lowpass or highpass.');

        end
    else

        filt = filt + 1;% Everything equals to one --> no damping
        map_dataf = bg;
        disp('You have set llp and hpp to zero! No filtering is applied')
    end
        

    map_dataf = map_dataf(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)));
%     filt = filt(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)));

end

