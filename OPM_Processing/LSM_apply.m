function Ims = LSM_apply(Ims,V,C,B)
% Apply LSM to imates
%
% Input:
% V = templates (y-x-num_templates)
% C = coefficient matrix of convolution of templates with filter
% B = Butterworth mask
%
%% Read input

[y_dim,x_dim,num_im] = size(Ims);
num_tpl = size(V,3);

%% FOR EACH IMAGE: Coefficient matrix and solve system

for im_ii = 1:num_im
        
    % ImPC contains: (Im * Vj) x B
    Cim = zeros(num_tpl,1,y_dim*x_dim);
    for ii = 1 : num_tpl
        tmp = convB(Ims(:,:,im_ii).*V(:,:,ii),B);
        Cim(ii,1,:) = tmp(:);
    end
    
    % Solve system
    coeffs = lin2solv(cat(2,C,Cim));
    
    % remove artifacts and store
    coeffs = reshape(coeffs.',y_dim,x_dim,num_tpl);
    Ims(:,:,im_ii) = Ims(:,:,im_ii) - sum(coeffs.*V,3);
    
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

function s = lin2solv(mat)
%	Solves an array of linear equations implementing the naive Gauss
%	algorithm e.g Ax = B : A ^ B => x
%	Assumes matrix = [A B], with A = cat(3,A1,...,An) and likewise for B.
%	Is intended for large measurement arrays i.e. no purely singular entries

%	SYNTAX : solution = lin2solv(matrix)

%	Tomer Fekete, 2006
[row,col,~] = size(mat);
for ii = 1 : col - 1
    mat(ii,:,:) = mat(ii,:,:)./mat(ii,ii*ones(1,col),:);
    for jj = ii+1 : row
        mat(jj,:,:) = mat(jj,:,:)-mat(ii,:,:).*mat(jj,ii*ones(1,col),:);
    end
end
for jj = row : -1 : 2
    mat(1:jj-1,jj,:) = mat(1:jj-1,jj,:).*mat(jj*ones(1,jj-1),end,:);
    mat(1:jj-1,end,:) = mat(1:jj-1,end,:) - mat(1:jj-1,jj,:);
end
s = squeeze(mat(:,end,:));

end
