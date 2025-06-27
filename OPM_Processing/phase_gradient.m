function [dtheta_dx, dtheta_dy] = phase_gradient(z, dx, dy)
% PHASE_GRADIENT Computes the spatial gradient of the phase of a 2D complex field.
%
% Inputs:
%   z   - complex 2D field (matrix)
%   dx  - spacing in x-direction
%   dy  - spacing in y-direction
%
% Outputs:
%   dtheta_dx - ∂θ/∂x, x-component of phase gradient
%   dtheta_dy - ∂θ/∂y, y-component of phase gradient

    % Extract real and imaginary parts
    psiR = real(z);
    psiI = imag(z);
    
    % Compute spatial derivatives
    [dpsiR_dx, dpsiR_dy] = gradient(psiR, dx, dy);
    [dpsiI_dx, dpsiI_dy] = gradient(psiI, dx, dy);
    
    % Compute |z|^2 = ψ_R^2 + ψ_I^2
    abs2 = psiR.^2 + psiI.^2;
    
    % Avoid division by zero
    abs2(abs2 == 0) = eps;

    % Gradient of the phase θ = arg(z)
    dtheta_dx = (psiR .* dpsiI_dx - psiI .* dpsiR_dx) ./ abs2;
    dtheta_dy = (psiR .* dpsiI_dy - psiI .* dpsiR_dy) ./ abs2;
end
