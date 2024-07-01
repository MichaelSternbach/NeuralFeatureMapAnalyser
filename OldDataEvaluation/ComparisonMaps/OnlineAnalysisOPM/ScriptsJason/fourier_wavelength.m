function [lambda, max_power] = fourier_wavelength(map)
%FOURIER_WAVELENGTH Calculate the average wavelength of a map.
%
% [lambda, max_power] = fourier_wavelength(map)
%
% Calculates the average wavelength, in pixels, of the given map. Calculated
% as the mean Fourier wavelength of the map averaged over all directions.
% The maximum Fourier power is also returned.
%
% Based on code originally written by Miguel A. Carreira-Perpinan.

% Input validation
if ~ismatrix(map)
    error('Input has too few or too many dimensions.');
end

G = size(map);
M = prod(G);

% Fourier spectrum
F = fft2(map);
F = real(F.*conj(F));

% Zero the DC component
F0 = F;
F0(1,1) = 0;
% Shift k = 0 to centre for display
F = fftshift(F);
F0 = fftshift(F0);
% Flatten it
f = reshape(F,M,1);
f0 = reshape(F0,M,1);

[X,Y] = ndgrid(1:G(1),1:G(2));

% Normalised square-root spectra
p = f/sum(f);					% Original
p0 = f0/sum(f0);				% Zero-DC

% Origin of the Fourier space
Fm = p'*[X(:) Y(:)];			% Mean of the Fourier space

% Fourier space Cartesian coordinates
X = X - Fm(1);
Y = Y - Fm(2);

% Power of dominant wavelength
max_power = max(f0/(sum(map(:)).^2));

X = X(:)/G(1);
Y = Y(:)/G(2);

% These are the moduli of the transformed (k1,k2), i.e., k = sqrt(k1^2+k2^2).
k = sqrt(sum([X Y].^2,2));

% Mean wave number
k(k == 0) = 1;                      % Where k=0, p0=0 so p0/k = 0
km = 1 / sum(p0./k);				% km is in cycles per pixel

% Mean wavelength
lambda = 1/km;                      % Pixels
