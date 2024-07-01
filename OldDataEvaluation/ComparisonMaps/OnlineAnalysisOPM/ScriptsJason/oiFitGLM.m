function [Z, thresh, x] = oiFitGLM(images, A, c, varargin),
%OIFITGLM Fit a General Linear Model (GLM).
%
%  [Z, THRESH] = OIFITGLM(IMAGES, A, C[, opt1, val1, ...]) fits a GLM
%  defined by A (the design matrix) to the supplied image sequence.
%
%  Available options are:
%   

% FIXME: add options as appropriate
%   'xxx'       - a vector specifying something...

% $Id: $

% % check input and output arguments
% error(nargchk(1, 4, nargin));



% default options...
frameRate = 5.0; % Hz
stimOnset = 1.0; % seconds

sigma = 2.0; % sdt dev of the Gaussian filter
p = 0.01; % significance threshold for the SPM{Z}'s

if nargin > 3,
  if isstruct(varargin{1}),
    opts = varargin{1};
  else,
    if ischar(varargin{1}),
      % parse the options strings
      opts = parseOpts({'frameRate',frameRate,'stimOnset',stimOnset,'sigma',sigma,'p',p}, varargin{:});
%       opts = parseOpts({'sigma',sigma,'p',p}, varargin{:});
    else,
      error('Too many input arguments, and/or no invalid options found. Type "help oiFitGLM" for usage information.');
    end
  end
     
  % process the options structure
  if isfield(opts,'frameRate'), ...
      frameRate = opts.frameRate; end
  if isfield(opts,'stimOnset'), ...
      stimOnset = opts.stimOnset; end
  if isfield(opts,'sigma'), ...
      sigma = opts.sigma; end
  if isfield(opts,'p'), ...
      p = opts.p; end
end

framePeriod = 1.0/frameRate;

% General Linear Model
[numConds,numTrials] = size(images);

[m,n,numFrames] = size(images{numConds,numTrials});

% spatial 'gaussian' filter
% h = fspecial('gaussian',5,sigma);
h = fspecial('gaussian',10,sigma);

for condId = 1:numConds,
  for cnt = 1:numTrials,
    for frame = 1:numFrames,
      images{condId,cnt}(:,:,frame) = imfilter(images{condId,cnt}(:,:,frame),h,'symmetric','conv');
    end
  end
end

% temporal filter/smoothing matrix
K = eye(numFrames); % no smoothing

% % construct design matrix
A = zeros([numFrames,4]); %6]);
% 
% % the hemodynamic response function (HRF)
A(stimOnset*frameRate,1) = 1;
tmp = conv(A(:,1), gampdf([0:0.2:10]',7.0,1));
A(:,1) = tmp(1:numFrames);
% 
t = [0:numFrames-1]'*framePeriod;
% 
% % MA017
% %% 0.55Hz sine and cosines - respiration
A(:,2) = cos(2*pi*0.55*t);
A(:,3) = sin(2*pi*0.55*t);
% 
% % 0.1Hz sine and cosines
% A(:,2) = cos(2*pi*0.1*t);
% A(:,3) = sin(2*pi*0.1*t);
% 
% % CatAF
% % 0.4Hz sine and cosines - respiration
% % A(:,4) = cos(2*pi*0.4*t);
% % A(:,5) = sin(2*pi*0.4*t);
% 
% % linear ramp
A(:,4) = t/max(t);

% smooth the design matrix
A = K*A;

% construct the "residual forming matrix"
R = eye(numFrames) - A*inv((A'*A))*A';

% estimate the no. of degrees of freedom
v = trace(R*K*K')^2/trace(R*K*K'*R*K*K');

c = [1, 0, 0, 0]; %, 0, 0]; % "select" the amplitude of the HRF basis function 

clear x
for condId = 1:numConds,
  for cnt = 1:numTrials,
    for ii = 1:m,
      for jj = 1:n,
        y = K*double(squeeze(images{condId,cnt}(ii,jj,:)));

        x{condId,cnt}(ii,jj,:) = inv(A'*A)*A'*y;
%         x{condId,cnt}(ii,jj,:) = A\y;
        
        % unbiased estimate of the variance of the noise
        s2{condId,cnt}(ii,jj) = (R*y)'*(R*y)/trace(R*K*K');

        % Z-scores (Z satisfies a t-distribution with v DOF; if v > 45 the
        % t-distribution is approximated by a Normal distribution)
        Z{condId,cnt}(ii,jj) = c*squeeze(x{condId,cnt}(ii,jj,:))/sqrt(c*s2{condId,cnt}(ii,jj)*inv(A'*A)*A'*K*K'*A*inv(A'*A)*c');
      end
    end
  end
end

% calculate the threshold for the SPM{Z}...
% u = 0:0.001:10;
% pz = (m*n)*(2*pi)^(-3/2)*2^(-1)*sigma^(-2)*u.*exp(-u.^2/2);
% thresh = u(max(find(pz > p)));
thresh = oiCalcThresh(p, m*n, sigma);
