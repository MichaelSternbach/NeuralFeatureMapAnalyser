function Ims = LSM(V,Ims,filt_w,templates)
%{
Implementation of the local similarity minimization (LSM) algorithm. 

Paper:

Removal of spatial biological artifacts in functional maps by local similarity minimization.
Fekete T1, Omer DB, Naaman S, Grinvald A.

J Neurosci Methods. 2009 Mar 30;178(1):31-9. doi: 10.1016/j.jneumeth.2008.11.020. Epub 2008 Dec 3.

%}
% If only one argument, display components
if nargin==1
   
    [y_dim,x_dim,~] = size(V);
    
    % calculate templates from singular value decomposition
    V = reshape(V,x_dim*y_dim,size(V,3));
    %[V,~,~] = svd(V,'econ');V = V(:,1:num_tpl); % <- slower
    [e_vec,e_val]=eig(V'*V);
    [e_val,ind]=sort(diag(e_val),'descend');
    e_vec=e_vec(:,ind);
    V=V*e_vec(:,1:16)*diag(e_val(1:16).^-0.5);
    V = reshape(V,y_dim,x_dim,16);
    
    % make figure
    figure
    for ii=1:1:16
        subplot(4,4,ii)
        imagesc(V(:,:,ii))
        title(ii)
        axis image off
    end

    Ims = [];
    return
end

%%
% wsiz : size of window to match template
% pcNum : number of templates
% templates : list of eigenimages to take from the blank

[y_dim,x_dim,num_im] = size(Ims);
if numel(templates)==1
    templates = 1:templates;
end
num_tpl = length(templates);

%% FOR ALL: Filter and Coefficient matrix

% creates a wsiz padded Butterworth mask
ord = 3;
s = filt_w./((1/0.95-1).^(1/(2*ord)));
[X,Y]= meshgrid(((0:x_dim-1)- round((x_dim-1)/2)),((0:y_dim-1)- round((y_dim-1)/2)));

B = 1./(1+(sqrt(X.^2+Y.^2)./s).^(2*ord));
B = fft2([zeros(filt_w,x_dim+2*filt_w);...
    zeros(y_dim,filt_w) B zeros(y_dim,filt_w);...
    zeros(filt_w,x_dim+2*filt_w)]);

% calculate templates from singular value decomposition
V = reshape(V,x_dim*y_dim,size(V,3));
%[V,~,~] = svd(V,'econ');V = V(:,1:num_tpl); % <- slower
[e_vec,e_val]=eig(V'*V);           
[e_val,ind]=sort(diag(e_val),'descend');   
e_vec=e_vec(:,ind);
V=V*e_vec(:,templates)*diag(e_val(templates).^-0.5);  
V = reshape(V,y_dim,x_dim,num_tpl);

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

%% FOR IMAGE: Coefficient matrix and solve system

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
