function C = corr_matrix(X, dr)
%CORR_MATRIX Calculate a shifted cross-correlation matrix.
%
% C = corr_matrix(X, dr)
%
% Calculates a shifted cross-correlation function for images across the
% third-dimension of the input matrix X. The input dr = [dr_x, dr_y] defines
% the shift. A value of dr = [0, 0] will perform the same function as the
% built-in function cov.
%
% Follows the definition of the cross-correlation function given in equation
% 4.56 of:
% Schiessl I (2001) "Blind source separation algorithms for the analysis of
% optical imaging experiments." PhD Thesis, Technische Universität Berlin.
% http://webdoc.sub.gwdg.de/ebook/diss/2003/tu-berlin/diss/2001/schiessl_ingo.pdf
%
% See also:
% COV
% ESD

if ndims(X) == 3
    [ny, nx, N] = size(X);
    n = nx*ny;
else
    error('X must be 3-dimensional.');
end

if all(dr == [0,0])
    C = cov(reshape(X, n, N));
else
    C = zeros(N, N);
    for i = 1:N
        Y = imtranslate(X(:,:,i), dr, 'FillValues', NaN);
        Q = sum(~isnan(Y(:)));
        for j = 1:N
            V = Y.*X(:,:,j);
            C(i,j) = nansum(V(:))/Q;
        end
    end
end