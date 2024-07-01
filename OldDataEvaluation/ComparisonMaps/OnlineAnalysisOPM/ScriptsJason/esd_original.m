
function [Y, W] = esd_original(X, dr)
%ESD Perform extended spatial decorrelation.
%
% [Y, W] = esd(X, dr)
%
% Performs extended spatial decorrelation, a form of blind-source separation,
% on the mixed 2D source images in matrix X (of size width*height*num_sources).
% The input dr = [dr_x, dy_y] is the 2D shift to use for the shifted
% cross-correlation matrix. The output Y is a matrix the same size as X
% containing the separated sources, and W is the demixing matrix.
%
% The coefficients of each separated source in the original data are given by
% the columns of inv(W).
%
% Follows the method described in section 4.5.1 of:
% Schiessl I (2001) "Blind source separation algorithms for the analysis of
% optical imaging experiments." PhD Thesis, Technische Universität Berlin.
% http://webdoc.sub.gwdg.de/ebook/diss/2003/tu-berlin/diss/2001/schiessl_ingo.pdf

% Handle input
if nargin < 2
    dr = [5, 5];
end

% Settings
height = size(X, 1);
width = size(X, 2);
N = size(X, 3);
px = width;
py = height;
n = px*py;


% Calculate zero-shift correlation matrix
C0 = corr_matrix(X, [0, 0]);
C0 = 0.5*(C0 + C0.');

% Derive the eigenvalues equation
[V, L] = eig(C0);
D = L^(-1/2)*(V.');

% Sphere the data
Y = reshape((D*reshape(X, n, N).').', [py, px, N]);

% Calculate shifted correlation matrix
Cdr = corr_matrix(Y, dr);
% Artifically symmetrise the matrix (as Cdr is often not perfectly symmetric
% for real world datasets)
Cdr = 0.5*(Cdr + Cdr.');

% Find the demixing matrix W
[V, ~] = eig(Cdr);
U = V.';
W = real(U*D);

% Returned the separated sources
Y = reshape((W*reshape(X, n, N).').', [py, px, N]);