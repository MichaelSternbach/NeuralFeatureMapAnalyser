function [average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient] = get_column_spacingManuel(map,ROI,pixel_per_mm,smallest_w,largest_w,w_step,interpolation_method)

if nargin<7 
	interpolation_method='polynomial';
end 

smallest_w = smallest_w * 1000; 
largest_w = largest_w * 1000; 
w_step = w_step * 1000; 
absolute_scale = 1000/pixel_per_mm; %% micrometers per pixel
%calculating the wavelength in the desired area - deleted the shrinkage
%information not sure what this was for in the first place
 disp('Estimate local column spacing of real part via wavelet analysis ...');
[average_w_re, local_w_re, new_roi_re,WavletCoefficient_re] = wavelength_estimator_data(real(map), ROI, smallest_w, largest_w, w_step,absolute_scale, interpolation_method);
disp('Estimate local column spacing of imaginary part via wavelet analysis ...');
[average_w_im, local_w_im, new_roi_im,WavletCoefficient_im] = wavelength_estimator_data(imag(map), ROI, smallest_w, largest_w, w_step,absolute_scale, interpolation_method);
average_spacing_mm = (average_w_re+average_w_im)/2/1000;

if nargout >3
   WavletCoefficient.X = (WavletCoefficient_re.X+WavletCoefficient_im.X)/2;
   WavletCoefficient.Y_mean = (WavletCoefficient_re.Y_mean+WavletCoefficient_im.Y_mean)/2;
   WavletCoefficient.XI = (WavletCoefficient_re.XI+WavletCoefficient_im.XI)/2;
   WavletCoefficient.YI_mean =(WavletCoefficient_re.YI_mean+WavletCoefficient_im.YI_mean)/2;
end

if nargout >1
    local_spacing_mm = (local_w_re+local_w_im)/2/1000;
    local_spacing_mm((new_roi_re.*new_roi_im) == 0) = 0;
end
if nargout >2
    newROI = new_roi_re.*new_roi_im;
end

disp('DONE.');

end

