function [f,phi,SN,sn] = GIF(f,SN_th)
% Implementation of the 'Generalized Indicator Functions' method by Yokoo
% et al.  see http://www.ncbi.nlm.nih.gov/pubmed/11707087 for Equations.
%
% GIF (phi) are a family of images that:
%   - maximize variations between conditions (signal) Eq.8
%   - minimize variations inside each condition (noise) Eq.11
% Projections of the data to a set of GIFs with signal/noise above a 
% threshold (SN_th) are returned.
%
% Since the optimization problem (Eq. 18) is composed of a linear combination
% of the original images, the snapshot method is used. The data is
% decomposed in singular values (Eq. 27) and an equivalent optimization
% problem is constructed (Eq. 38). The solutions 'beta' are then transformed
% back to image space.
%
% Input:
% f = 4d array of grayscale images with dimensions
%       1:pixel-y  2:pixel-x  3:condition (m from M)  4:block (t from T)
% SN_th = signal/noise threshold for GIFs
%
% Output:
% f = Projections of GIF on individual images
% phi = used GIFs
% SN = signal/noise of the used GIFs

%% INPUT

% determine the signal/noise threshold
if nargin<2
    SN_th = 3;
end

% get dimensions of data [M = conditions, T = num_blocks]
[pix_y,pix_x,M,T] = size(f);

% remove mean over blocks and conditions
f_ref = mean(mean(f,4),3);
f = bsxfun(@minus,f,f_ref);
 
% reshape image in 2D 1:pixel 2:image (M*T, ordered in conditions)
f = reshape(permute(f,[1 2 4 3]),pix_y*pix_x,T*M);

%% GIF METHOD

% === (Eq. 27)
% decompose the images into singular values (PCA)
[a_mt,sigma]=eig(f'*f);
[sigma,idx] = sort(diag(sigma),'descend');
a_mt=a_mt(:,idx);
% psi = f*a*diag(sigma.^-0.5); % <- not stored to save memory

% === (Eq. 1,2)
% do averages for the next calculations
a_m_t = permute(reshape(a_mt',M*T,T,M),[3 2 1]);
a_m = squeeze(mean(a_m_t,2));
a = mean(a_m,1);

% === (Eq. 40):7
% get Signal Covariance Matrix 
Cs = zeros(M*T,M*T);
for m=1:M
    tmp = bsxfun(@minus,a_m(m,:),a);
    Cs = Cs + tmp'*tmp;
end
Cs = (1/(M-1))*Cs;

% === (Eq. 41):10
% get Noise Covariance Matrix
Cn = zeros(M*T,M*T);
for m =1:M
    tmp = bsxfun(@minus,squeeze(a_m_t(m,:,:)),a_m(m,:));
    Cn  = Cn + tmp'*tmp;
end
Cn = (1/(M*T-M))*Cn;

% === (Eq. 42)
% diagonal matrix of eigenvalues
D = diag(sigma.^(0.5));

% === (Eq. 38):18
% solve eigenvalue problem 
[beta,sn] = eig(D*Cs*D - SN_th*(D*Cn*D));

% === (Eq. 23)
% identify solutions above signal/noise level
beta = beta(:,diag(sn)>0);
if ~isdeployed
    %disp([num2str(size(beta,2)),' components are above SN threshold']);
end
if isempty(beta)
    f =[];phi=[];
    return;
end

% == (Eq. 3)
% response amplitude functions in svd space (same as in image space)
rho = a_mt*D*beta;

% == (Eq. 37)
% get the GIFs in image space [NOTE: psi = f*a*diag(sigma.^-0.5)]
phi   = (f*a_mt*diag(sigma.^-0.5))*beta; 

% == (Eq. 3)
% use same response amplitude functions on phi
% project the GIF to the data
f = rho*phi';

%% OUTPUT

% == reshape output and add mean
f = permute(reshape(f',pix_y,pix_x,T,M),[1 2 4 3]);
f = bsxfun(@plus,f,f_ref);

% average
%f = mean(f,4);

if nargout>1
    % == reshape GIF
    phi = reshape(phi,pix_y,pix_x,size(beta,2));
end

if nargout>2
    % == (Eq. 20)
    % signal/noise of selected GIF
    rho = reshape(rho,T,M,size(beta,2));
    SN = squeeze(var(mean(rho,1),[],2)./mean(var(rho,[],1),2));
end

end