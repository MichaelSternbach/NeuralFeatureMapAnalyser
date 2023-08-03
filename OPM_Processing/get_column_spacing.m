function [average_spacing_mm,local_spacing_mm] = get_column_spacing(map,ROI,pixel_per_mm,smallest_w,largest_w,w_step,N_orientations)
% Get local spacing in the map
%
% INPUT:
%
% z : filtered orientation map
% ROI : region of interest
% pixels_per_mm : measure of the recoding
% smallest_w:w_step:largest_w : list of wavelet wavelenghts to test
% N_orientations : 0 if circular symmetric, N for sampling waveforms
%
% OUTPUT:
%
% average_spacing_mm : average spacing over new ROI
% local_spacing_mm : column spacing for each pixel in mm
% ROI_new : ROI where the column spacing could be estimated

%% read input

% N_orientations = 0 if circular symmetric, N for sampling waveforms
if nargin==6 
    N_orientations = 16;
end
k_base = 7; % roughly corresponds to no. waves within Gaussian
aniso = 1; % anisotropy parameter, for bandedness estimation, this parameter can be set to values larger than one to achieve anistropic wavelets
% (see Keil et al. PNAS, 2010 for more details)
wavelengths = smallest_w:w_step:largest_w;

%% Prepare padding for column spacing

% Determine size of padding
map_size=size(map);
max_size_needed  = max(map_size);
if max_size_needed < (10*largest_w*pixel_per_mm)
    max_size_needed = 10*largest_w*pixel_per_mm;
end
padd_size =round(2.^(ceil(log(max_size_needed)/log(2))+1));

% get indices of where the original data is stored
roi_padd = false(padd_size);
roi_padd(round((padd_size-map_size(1))/2+1):round((padd_size-map_size(1))/2+map_size(1)),round((padd_size-map_size(2))/2+1):round((padd_size-map_size(2))/2+map_size(2)))=ROI;

% zero padding and fourier transform or real and imaginary part
map_padd_re=zeros(padd_size);
map_padd_re(round((padd_size-map_size(1))/2+1):round((padd_size-map_size(1))/2+map_size(1)),round((padd_size-map_size(2))/2+1):round((padd_size-map_size(2))/2+map_size(2)))...
    =real(map);
map_padd_re = fft2(map_padd_re);

map_padd_im=zeros(padd_size);
map_padd_im(round((padd_size-map_size(1))/2+1):round((padd_size-map_size(1))/2+map_size(1)),round((padd_size-map_size(2))/2+1):round((padd_size-map_size(2))/2+map_size(2)))...
    =imag(map);
map_padd_im = fft2(map_padd_im);

% prepare filter size using max
length_scaled = smallest_w*pixel_per_mm * k_base/(2*pi);
N_pix = 2*ceil(8*length_scaled*aniso/2);
[X,Y] = meshgrid( (1:N_pix) - ceil(N_pix/2) );

filter_padd=zeros(padd_size);
filter_padd((padd_size-N_pix)/2+1:(padd_size-N_pix)/2+N_pix,(padd_size-N_pix)/2+1:(padd_size-N_pix)/2+N_pix)=1;
filter_padd_idx = find(filter_padd == 1);

%% Get column spacing

disp('CALCULATING TYPICAL SCALE.');
disp('   Doing convolutions with wavelets...');

% variables for storing results of wavelets (results in columns)
values_re = zeros(sum(ROI(:)),length(wavelengths));
values_im = zeros(sum(ROI(:)),length(wavelengths));

% loop over sizes
for wave_ii = 1:length(wavelengths)
    
    % Coordinate system
    lambda_pixel = wavelengths(wave_ii)*pixel_per_mm;
    length_scaled = lambda_pixel * k_base/(2*pi);
    
    if N_orientations == 0
        % ====== Rot Symmetric <  PLUG HERE OTHER WAVELET FILTERS!!!
        
        % make Mexican hat wavelet
        filter = 1/(pi*length_scaled^4)*(1-(X.^2+Y.^2)/(2*length_scaled^2)).*exp(-(X.^2 + (Y/aniso).^2)/(2*length_scaled^2));
        filter = filter - mean(filter(:));
        filter_padd(filter_padd_idx)=filter(:);
        
        % Convolute in Fourier domain and store
        tmp = abs(ifft2(map_padd_re.* fft2(fftshift(filter_padd))));
        values_re(:,wave_ii) = tmp(roi_padd);
        
        tmp = abs(ifft2(map_padd_im.* fft2(fftshift(filter_padd))));
        values_im(:,wave_ii) = tmp(roi_padd);
        
    else
        % ====== Gabor patches
        
        % Loop over filter orientations
        for theta = (0:N_orientations-1)/N_orientations*pi
            
            % Rotate coordinate system clockwise
            X_rot = cos(theta)*X + sin(theta).*Y;
            Y_rot = -sin(theta)*X + cos(theta).*Y;
            
            % Make Gabor wavelet and padd
            filter = ( cos( 2*pi*X_rot/lambda_pixel ) + 1i*sin( 2*pi*X_rot/lambda_pixel ) ) ...
                .*(1/length_scaled*exp(-(X_rot.^2 + (Y_rot/aniso).^2)/(2*length_scaled^2)));
            filter = filter - mean(filter(:));
            filter_padd(filter_padd_idx)=filter(:);
            
            % Convolve in Fourier domain and store
            tmp = abs(ifft2(map_padd_re.* fft2(fftshift(filter_padd))));
            values_re(:,wave_ii) = values_re(:,wave_ii) + tmp(roi_padd)/N_orientations;
            
            tmp = abs(ifft2(map_padd_im.* fft2(fftshift(filter_padd))));
            values_im(:,wave_ii) = values_im(:,wave_ii) + tmp(roi_padd)/N_orientations;
            
        end
    end
end

% empty space before fitting
clearvars map_re map_im filter_padd

%% interpolate real and imaginary part separately and find local maxima

disp('   Finding column size for each pixel...');
% interpolate values
k_interp = 0.005;
XI = smallest_w:w_step*k_interp:largest_w;

% valid range
XI_range = [find(XI>(smallest_w+w_step),1) find(XI<(largest_w-w_step),1,'last')];

max_peak_re = zeros(sum(ROI(:)),2);
max_peak_im = zeros(sum(ROI(:)),2);

for pix_ii=1:sum(ROI(:))
    % real part
    values_re_ii = interp1(wavelengths,values_re(pix_ii,:),XI, 'spline');
    idx = 1+find(values_re_ii(2:end-1)>=values_re_ii(1:end-2) & values_re_ii(2:end-1)>=values_re_ii(3:end));
    idx(idx<XI_range(1) | idx>XI_range(2)) = [];
    if ~isempty(idx)
        [val,tmp] = max(values_re_ii(idx));
        max_peak_re(pix_ii,:) = [XI(idx(tmp)),val];
    end
    % imaginary part
    values_im_ii = interp1(wavelengths,values_im(pix_ii,:),XI, 'spline');
    idx = 1+find(values_im_ii(2:end-1)>=values_im_ii(1:end-2) & values_im_ii(2:end-1)>=values_im_ii(3:end));
    idx(idx<XI_range(1) | idx>XI_range(2)) = [];
    if ~isempty(idx)
        [val,tmp] = max(values_im_ii(idx));
        max_peak_im(pix_ii,:) = [XI(idx(tmp)),val];
    end
end

% do weighted average
max_peak = (max_peak_re(:,2).*max_peak_re(:,1) + max_peak_im(:,2).*max_peak_im(:,1))./(max_peak_re(:,2)+max_peak_im(:,2));
max_peak(isnan(max_peak)) = 0;

% Convert to 2D
local_spacing_mm = zeros(size(ROI));
local_spacing_mm(ROI) = max_peak;
average_spacing_mm = mean(local_spacing_mm(local_spacing_mm>0));

disp('DONE.');

end

