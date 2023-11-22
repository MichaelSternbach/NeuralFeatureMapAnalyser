function [average_w, local_w, new_roi,WavletCoefficient] = wavelength_estimator_data(map, roi_matrix, smallest_w, largest_w, w_step,absolute_scale, interpolation_method)
%
% FUNCTION CALLED BY ANALYZE_SINGLE_MAP.M
%
% USAGE : 
%     wavelength_estimator(map_data, roi_matrix,small_scale, large_scale, scale_step, absolute_scale)
%     This function estimates the local wavelength of an OPM
%     The map is zero padded for the FFT, and after computing the filtering
%     resized to the original size
%     if wavelet coefficients don't show a clear maximum
%     within the scales scanned by the routine, such parts of the data a
%     ormitted from wavelength estimation
%
%
% INPUT PARAMETERS
% map_data      ... preprocessed complex-valued OPM data 
% roi_matrix    ... region of interest, matrix with 1 inside the ROI, and
%                    zero otherwise
%
% smallest_w : smallest wavelength possible (in micrometers)
% largest_w : largest wavelength possible (in micrometers)
% w_step :  step (in micrometers) Kaschube et al. Science 2010, 50mu m
%
% absolute scale ... micrometers per pixel for the given data matrix
% interpolation_method ... string, choose either 'polynomial',
%                          'spline' or 'none'
%                           default is polynomial (see Keil et al. PNAS, 2010)
%
% OUTPUT PARAMETERS
% average_w         ... scalar value, average local wavelength in pixels
% local_w           ... matrix, same size as map_data with local wavelength
%                       according to the wavelet estimate
% new_roi           ... new ROI matrix, only containing ones where local
%                       column spacing could be estimated
%
%  Copyright (C) 2014, Max-Planck-Institut for Dynamics and Self-organization, The  
%  Nonlinear Dynamics Group. This software may be used, copied, or 
%  redistributed as long as it is not sold and this copyright notice is 
%  reproduced on each copy made. This routine is provided as is without 
%  any express or implied warranties whatsoever.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    scrsz = get(0,'ScreenSize');

    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    
    if nargin < 7
        interpolation_method  = 'spline';
    end
    
    
    k_base = 7;% Wavelet parameters, k_base roughly corresponds to no. waves within Gaussian
    aniso = 1; % Anisotropy parameter, set to one, 
               % for bandedness estimation, this parameter can be set to values larger than one to achieve anistropic wavelets
               % (see Keil et al. PNAS, 2010 for more details) 
    
    
    a=size(map);

    max_size_needed  = max(size(map));
    if max_size_needed < (10 * largest_w/absolute_scale)
        max_size_needed = (10 * largest_w/absolute_scale);
    end
    
    
    %Find next power of two --> size of the padded array
    ext =round(2.^(ceil(log(max_size_needed)/log(2))+1));
    
    % put map on square array of size 2^N (zero padding)
    map_data_padd=zeros(ext);     

    map_data_padd(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)))=map;
    
        
    % Fourier Transform the image in order to convolute it later on
    map_data_padd = fft2(fftshift(map_data_padd));
    
    N = 16; % number of different orientations of wavelet patches are used
    
    components = zeros(N, ext,ext); % Field for the components of the individually rotated filters
    orient_av_mod = zeros(ext,ext,length(smallest_w:w_step:largest_w));% Field for the average value for one scale    


    
    index = 0;    
    % Loop over scales (wavelengths)
      
    for wavelength = smallest_w:w_step:largest_w
        index = index+1;

        % Loop over filter orientations
        for j = 0:1:N-1

            theta = pi/N*j;
            Lambda_pixel = wavelength/absolute_scale;
            k_pixel = 2*pi/Lambda_pixel;
            l = k_base/k_pixel;
            % Generate a filter with 8 sigma width and even size matrix
            filter = generate_wavelet_filter(round(8*l*aniso) + mod(round(8*l*aniso),2),wavelength,k_base,aniso,theta,absolute_scale);
            
            %%% Padd filter for FFT
            filter_padd=zeros(ext);
            b = size(filter);
            filter_padd((ext-b(1))/2+1:(ext-b(1))/2+b(1),(ext-b(2))/2+1:(ext-b(2))/2+b(2))=filter;            
             
            %%%%% Filtern in Fourier Domain
            filter = fft2(fftshift(filter_padd));
            filt_im = abs(fftshift(ifft2(map_data_padd.* filter)));
                        
            components(j+1,:,:) = filt_im;
            
        end % For loop over orientations
        
        orient_av_mod(:,:, index) = 1/pi/N*sum(abs(components),1);
    end% For loop over scales

    % Resize the field (reverse zero padding)
    orient_av_mod = orient_av_mod(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)),:);
    disp(['Computation of wavelet coefficients finished.']);
 
 

    %%%% Fit polynomials/splines to the wavelet coefficients
    k_interp = 0.005;
    X = smallest_w:w_step:largest_w;
    XI = smallest_w:w_step*k_interp:largest_w;
    
    Y = zeros(length(X),1);
    local_w = zeros(size(orient_av_mod,1),size(orient_av_mod,2));
    

    if strcmpi(interpolation_method, 'spline') || strcmpi(interpolation_method, 'splines')
        disp('Now interpolating wavelet coefficients with cubic splines ...');
    elseif strcmpi(interpolation_method, 'polynomial') ||  strcmpi(interpolation_method, 'polynomials')
        disp('Now interpolating wavelet coefficients with polynomial ...');
    elseif strcmpi(interpolation_method, 'none')
        disp('Using wavelet coefficients without interpolation...');
    end

    counter = 0;
    
    
    figure(3);
    set(gcf,'Position',[2*scrsz(3)/3 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3]); %[left, bottom, width, height]:
    
    
    Y_mean = zeros(length(X),1);
    YI_mean = zeros(1,length(XI));
    for i = 1:size(orient_av_mod,1)
        for j = 1:size(orient_av_mod,2)
            
            if roi_matrix(i,j) == 1 %% Do that only for the roi parts
                counter = counter + 1;
                Y(:,1) = orient_av_mod(i,j,:);

                % Choose spline or polynomial interpolation
                if strcmpi(interpolation_method, 'spline')
                    YI = interp1(X, Y, XI, 'splines' );
                elseif strcmpi(interpolation_method, 'polynomial')
                    P = polyfit(X,Y',6);
                    YI = P(1)*XI.^6 + P(2)*XI.^5 + P(3)*XI.^4 + P(4)*XI.^3 + P(5)*XI.^2 + P(6)*XI.^1 + P(7);
                elseif strcmpi(interpolation_method, 'none')
                    YI = Y; %% no interpolation
                end
        
                if counter < 100
                    figure(3);
                    plot(X, Y,XI,YI);
                    if strcmpi(interpolation_method, 'spline')
                        title('Wavelet coefficients and spline interpolation');
                    elseif strcmpi(interpolation_method, 'polynomial')
                        title('Wavelet coefficients and polynomial interpolation');
                    elseif strcmpi(interpolation_method, 'none')
                        title('Wavelet coefficients');
                    end
                    
                    title('Wavelet coefficients');
                    xlabel('Scale [mu m]');
                    ylabel('Wavelet Coeff');
                    set(gca, 'fontsize', 16);
                    drawnow
                    
                end
                
                Y_mean = Y_mean + Y;
                YI_mean = YI_mean + YI;

                index = find(YI == max(YI));
                if length(index) > 1 
                    index = index(1);
                end
                local_w(i,j) = XI(index);              
            end
        end
    end
    
    local_w = local_w.*roi_matrix;
    
    Y_mean = Y_mean./counter;
    YI_mean = YI_mean./counter;
    %% plot mean wavelet coefficient
    figure(4);
    plot(X, Y_mean,XI,YI_mean);
    title('mean Wavelet coefficients');
    xlabel('Scale [mu m]');
    ylabel('Wavelet Coeff');
    set(gca, 'fontsize', 16);
    
    if nargout >3
       WavletCoefficient.X = X;
       WavletCoefficient.Y_mean = Y_mean;
       WavletCoefficient.XI = XI;
       WavletCoefficient.YI_mean =YI_mean;
    end

    % Now check if we are too close to the boundaries of the scales
    new_roi = roi_matrix;
    q = find((local_w < (smallest_w + 2*w_step)) & local_w ~= 0);
    if ~isempty(find(q,1))
        disp('        Some regions have been deleted from the region of interest due to too small wavelength!');
        local_w(q) = 0;
        new_roi(q) = 0;
    end
    q = find(local_w > (largest_w - 2*w_step));
    if ~isempty(find(q,1))
        disp('        Some regions have been deleted from the region of interest due to too large wavelength!');
        local_w(q) = 0;
        new_roi(q) = 0;
    end
    average_w = mean(mean(local_w(new_roi == 1)));
%     
    disp('... wavelet analysis done.');
        
end %End of function WAVELET_ESTIMATOR_DATA.M




%%%%%%%%%%%%%%%%%%%%%%%%%% LIST OF USED FUNCTIONS %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wavelet = generate_wavelet_filter(N,wavelength,k,aniso,theta,absolute_scale)
%
% function wavelet = generate_wavelet_filter(N,wavelength,k,aniso,theta,absolute_scale)
%   USAGE:
%   INPUT PARAMETERS
%   N           ... size of the filter (is a square matrix)
%   wavelength  ... wavenumber in the filter (should be integer) attention to periodicity!!!!!!!!!!
%   aniso       ... anisotropy between x and y direction
%   theta       ... angle of the modulating wave
%   OUTPUT PARAMETERS: a complex-valued Morlet wavelet function


    Lambda_pixel = wavelength/absolute_scale;
    k_pixel = 2*pi/Lambda_pixel;
    [X Y] = meshgrid(1:N,1:N);
    X=X-ceil(N/2);
    Y=Y-ceil(N/2);

    wave = cos(cos(theta)*X*k_pixel + sin(theta)*Y*k_pixel) + 1i*sin(cos(theta)*X*k_pixel + sin(theta)*Y*k_pixel);

    l = k/k_pixel;
    
    %%% Rotated coordinate system
    X_rot = cos(theta)*X + sin(theta).*Y;
    Y_rot = -sin(theta)*X + cos(theta).*Y;

    gaussian = 1/l*exp(-1/2/l/l*(X_rot.*X_rot + 1/aniso/aniso*Y_rot.*Y_rot)); % Gives an anisotropic gaussian

    % Now multiply the wave pointwise with the gaussian
    wavelet = gaussian.*wave;    
    wavelet = wavelet - mean(mean(wavelet));
        
end




