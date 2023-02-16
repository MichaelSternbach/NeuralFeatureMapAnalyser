function [V,C,B] = LSM_get_templates(V,filt_w,templates)
% Get the templates V for LSM and calculate coefficient matrix C and Butterworth mask B
%
% Input:
% V = blanks 3D (y-x-blocks)
% filt_W = width of the filter
% templates = list or number of templates to take
%
%%

% If only one argument, display components
if nargin==1
   
    [y_dim,x_dim,blocks] = size(V);
    
    templates = min([16 blocks]);
    
    % calculate templates from singular value decomposition
    V = reshape(V,x_dim*y_dim,blocks);
    %[V,~,~] = svd(V,'econ');V = V(:,1:num_tpl); % <- slower
    [e_vec,e_val]=eig(V'*V);
    [e_val,ind]=sort(diag(e_val),'descend');
    e_vec=e_vec(:,ind);
    V=V*e_vec(:,1:templates)*diag(e_val(1:templates).^-0.5);
    V = reshape(V,y_dim,x_dim,templates);
    
    % make figure
    figure
    for ii=1:1:templates
        subplot(4,4,ii)
        imagesc(V(:,:,ii))
        title(ii)
        axis image off
    end

    V = [];
    C = [];
    B = [];
    return
end

%% READ INPUT

[y_dim,x_dim,num_blanks] = size(V);
if numel(templates)==1
    templates = 1:templates;
end
num_tpl = length(templates);

%% Get Templates, Filter and Coefficient matrix

% calculate templates from singular value decomposition
V = reshape(V,x_dim*y_dim,num_blanks);
%[V,~,~] = svd(V,'econ');V = V(:,1:num_tpl); % <- slower
[e_vec,e_val]=eig(V'*V);           
[e_val,ind]=sort(diag(e_val),'descend');   
e_vec=e_vec(:,ind);
V=V*e_vec(:,templates)*diag(e_val(templates).^-0.5);  
V = reshape(V,y_dim,x_dim,num_tpl);

% creates a filt_w padded Butterworth mask
ord = 3;
s = filt_w./((1/0.95-1).^(1/(2*ord)));
[X,Y]= meshgrid(((0:x_dim-1)- round((x_dim-1)/2)),((0:y_dim-1)- round((y_dim-1)/2)));

B = 1./(1+(sqrt(X.^2+Y.^2)./s).^(2*ord));
B = fft2([zeros(filt_w,x_dim+2*filt_w);...
    zeros(y_dim,filt_w) B zeros(y_dim,filt_w);...
    zeros(filt_w,x_dim+2*filt_w)]);

% Get the coefficient matrix
% C contains: (Vi * Vj) x B
C = zeros(num_tpl,num_tpl,y_dim*x_dim);
for ii = 1 : num_tpl
    tmp = convB(prod(V(:,:,[ii ii]),3),B);
    C(ii,ii,:) = tmp(:);
    for jj = ii + 1 : num_tpl;
        tmp = convB(prod(V(:,:,[ii jj]),3),B);
        C(ii,jj,:) = tmp(:);
        C(jj,ii,:) = tmp(:);
    end
end

end

function Im = convB(Im,filter_ft_pad)
% convolve matrix with given padded and fourier transformed filter

[y_dim,x_dim] = size(Im);
filt_w = (size(filter_ft_pad,1)-y_dim)/2;

% pad image with +/-filt_w
Im = [zeros(filt_w,x_dim+2*filt_w);...
    zeros(y_dim,filt_w) Im zeros(y_dim,filt_w);...
    zeros(filt_w,x_dim+2*filt_w)];

% convolve images with filter
Im = real(fftshift(ifft2(fft2(Im).*filter_ft_pad)));

% remove pad
Im = Im(filt_w+1:(end-filt_w),filt_w+1:(end-filt_w));

end
