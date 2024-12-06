function profile = calculate_power_profile_3d(map, profile_scale_pixels)
    % This function makes a radial average of the power spectrum of a 3D map
    % to get a spectrum profile.

    % initiate variable
    profile = zeros(length(profile_scale_pixels), 1);

    % calculate padding (for 3D)
    pad_size_pixels = max(size(map));%2^nextpow2(max(size(map)));  % Pad to the next power of two for FFT performance

    % calculate 3D power spectrum
    disp('calculate 3D power spectrum')
    PS = abs(fftshift(fftn(map, [pad_size_pixels, pad_size_pixels, pad_size_pixels]))).^2;

    % create function (3D)
    [kx, ky, kz] = ndgrid(-pad_size_pixels/2:pad_size_pixels/2-1, ...
                          -pad_size_pixels/2:pad_size_pixels/2-1, ...
                          -pad_size_pixels/2:pad_size_pixels/2-1);
    PowerFun = scatteredInterpolant(kx(:), ky(:), kz(:), PS(:));

    % calculate radial profile (in 3D)
    disp('calculate radial profile ')
    N_Angles = 60;%101;
    phi = linspace(0, 2*pi, N_Angles);  % Azimuthal angle
    phi(end) = [];
    theta = linspace(0, pi, N_Angles);  % Polar angle

    for ind = 1:length(profile_scale_pixels)
        r = pad_size_pixels / profile_scale_pixels(ind);  % Compute radius

        % Create spherical coordinates for radial averaging
        [PHI, THETA] = meshgrid(phi, theta);
        X = r .* sin(THETA) .* cos(PHI);
        Y = r .* sin(THETA) .* sin(PHI);
        Z = r .* cos(THETA);

        % Evaluate the power at the spherical coordinates and average
        profile(ind) = mean(PowerFun(X(:), Y(:), Z(:)));
    end

end
