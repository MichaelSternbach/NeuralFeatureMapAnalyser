function thresh = oiCalcThresh(p, s, sigma),
%OICALCTHRESH Calculate significance threshold for SPM{Z}'s.
%
%  THRESH = OICALCTHRESH(P, S, SIGMA) calculates the threshold at the
%  level of significance requested (0 < P < 1.0), for SPM{Z}'s of size
%  S = M*N pixels, assuming spatial Gaussian smoothing with standard
%  deviation SIGMA.
%
%  Available options are:
%    none

% $Id: $

% check input and output arguments
error(nargchk(3, 3, nargin));

u = 0:0.001:20;

pz = s*(2*pi)^(-3/2)*2^(-1)*sigma^(-2)*u.*exp(-u.^2/2);

thresh = u(max(find(pz > p)));