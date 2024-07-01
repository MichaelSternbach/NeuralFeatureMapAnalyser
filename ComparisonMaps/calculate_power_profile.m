function profile = calculate_power_profile(map,profile_scale_pixels)
% This function makes a radial average of the power spectrum of a map to
% get a spectrum profile.
%%

% initiate variable
profile = zeros(length(profile_scale_pixels),1);

% calculate padding
pad_size_pixels = 2^nextpow2(max(size(map)));

% calculate power spectrum
PS = abs(fftshift(fft2(map,pad_size_pixels,pad_size_pixels))).^2;
PS = PS/mean(PS(:));

% create function
[kx,ky] = meshgrid(-pad_size_pixels/2:pad_size_pixels/2-1,-pad_size_pixels/2:pad_size_pixels/2-1);
PowerFun = scatteredInterpolant(kx(:),ky(:),PS(:));

% calculate radial profile
theta=linspace(0,2*pi,101);theta(end) = [];
for ind = 1:length(profile_scale_pixels)
    r = pad_size_pixels/profile_scale_pixels(ind);
    profile(ind) = mean(PowerFun(bsxfun(@times,r,cos(theta)),bsxfun(@times,r,sin(theta))));
end

end