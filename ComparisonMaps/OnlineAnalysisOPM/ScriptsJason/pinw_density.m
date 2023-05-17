
function [density_per_pixel, density_per_wl, density_per_mm] = pinw_density(op, mask)
%PINW_DENSITY Calculate the pinwheel density of an orientation preference map.
%
% [density_per_pixel, density_per_wl] = pinw_density(op, mask)
%
% Calculates the average number of pinwheels either per pixel or per square
% map wavelength. The average map wavelength is computed using the map's
% Fourier spectrum. A binary mask can be given to restrict the analysis.

% Input defaults
if nargin < 2
    mask = true(size(op));
end

fp= exp(1i*op/256*2*pi);

% Calculate OP wavelength
op_wl = mean([fourier_wavelength(real(fp)), fourier_wavelength(imag(fp))]);

% Mask
op(~mask) = nan;

% Find pinwheels in OP map
pinw = locate_pinwheels(op);
npinw = size(pinw,1);

% Calculate per pixel
density_per_pixel = npinw/sum(mask(:));

% Calculate per wavelength
density_per_wl = density_per_pixel*op_wl^2;
