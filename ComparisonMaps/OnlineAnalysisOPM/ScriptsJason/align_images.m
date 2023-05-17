function [aligned] = align_images(images, ref_image, shape, method)
%ALIGN_IMAGES Aligns a collection of images to a reference image.
%
% [aligned, offsets, all_ccs] = align_images(images, ref_frame, shape, method)
%
% Aligns a given collection of images to a reference image using spatial
% translation by maximising linear correlation. Images should be a 3-dimensional
% stack of images to align, and ref_image should be an image of the same size,
% with which each image in the stack will be aligned. If ref_image is not given,
% the images will be aligned to the first in the stack.
% Shape can be one of:
%  'full'  - (default) returns the full alignment with NaNs in the missing
%            pixels.
%  'same'  - returns the alignment the same size as the reference image.
%  'valid' - returns the smallest region overlapping all aligned images.
% To specify shpae but not a refernce image, set ref_image to [].
% Method can be one of:
%  'brute' - (default) evaluate all correlations across a large range of shifts
%            and pick the maximum.
%  'optim' - use lattice_gradient_descent to find maximum correlation (much
%            quicker but may fail).
%
% The output aligned contains the aligned images, sized according to the shape
% parameter. Offsets contains a list of the translations used to align each
% image. All_ccs contains all correlation coefficients calculated, if the
% brute-force method was used.

% Input defaults and validation
  
if nargin == 0 || nargin > 4
    error('Requires between one and four inputs.');
end
if ndims(images) ~= 3
    error('Image matrix must be 3-dimensional.');
end
if nargin < 2 || isempty(ref_image)
    ref_image = images(:,:,1);
end
if nargin < 3
    shape = 'full';
end
if nargin < 4
    method = 'brute';
end
if ~strcmp(method, 'brute') && nargout == 3
    error('All correlation coefficients only available for method ''brute''.');
end

[height, width, nimages] = size(images);

if strcmp(method, 'brute')
    xs = -16:16;
    ys = -16:16;
end

if height ~= size(ref_image, 1) || width ~= size(ref_image, 2)
    error('Reference frame must be the same size as the data.');
end

offsets = nan(nimages,2);
if strcmp(method, 'brute')
    all_ccs = nan(nimages,length(xs),length(ys));
    for j = 1:nimages
        ccs = nan(length(xs), length(ys));
        for xi = 1:length(xs)
            for yi = 1:length(ys)
                ccs(xi,yi) = correlation(ref_image, images(:,:,j), [xs(xi), ys(yi)]);
            end
        end
        all_ccs(j,:,:) = ccs;
        [~, I] = max(max(ccs));
        [~, J] = max(max(ccs.'));
        offsets(j,:) = [xs(J), ys(I)];
    end
elseif strcmp(method, 'optim')
    for j = 1:nimages
        x = lattice_gradient_descent(@(x)(-correlation(ref_image,images(:,:,j),x)), [0 0], 100);
        offsets(j,:) = x;
    end
else
    error('Method must be ''brute'' or ''optim''.')
end

max_offset = max(abs(offsets(:)));
margin = 2*max_offset;
if mod(margin, 2) == 0
    margin = margin + 1;
end
mid = floor(margin/2);
temp = nan(height+margin, width+margin, nimages);

for j = 1:nimages
    offx = offsets(j,1);
    offy = offsets(j,2);
    temp(1+mid+offy:mid+offy+height,1+mid+offx:mid+offx+width,j) = images(:,:,j);
end

if strcmp(shape, 'full')
    aligned = temp;
elseif strcmp(shape, 'same')
    aligned = temp(1+mid:mid+height,1+mid:mid+width,:);
elseif strcmp(shape, 'valid')
    s = regionprops(sum(~isnan(temp),3)==nimages, 'BoundingBox');
    bb = s(1).BoundingBox;
    aligned = nan(bb(4), bb(3), nimages);
    bb(3) = bb(3)-1;
    bb(4) = bb(4)-1;
    for j = 1:nimages
        aligned(:,:,j) = imcrop(temp(:,:,j), bb);
    end
else
    error('Shape must be ''full'', ''same'', or ''valid''.')
end


function cc = correlation(a, b, s)
% Calculates the linear correlation between images A and B after the latter is
% translated by S = [X, Y] pixels.

im0 = a;
im1 = b;
x = s(1);
y = s(2);
if x < 0
    im1 = im1(:,1+abs(x):end);
    im0 = im0(:,1:end-abs(x));
elseif x > 0
    im1 = im1(:,1:end-x);
    im0 = im0(:,1+x:end);
end
if y < 0
    im1 = im1(1+abs(y):end,:);
    im0 = im0(1:end-abs(y),:);
elseif y > 0
    im1 = im1(1:end-y,:);
    im0 = im0(1+y:end,:);
end
n = numel(im0);
im0 = im0(:);
im1 = im1(:);
im0 = im0 - sum(im0)/n;
im1 = im1 - sum(im1)/n;
cc = sum((im0./sqrt(sum(im0.^2,1))).*(im1./sqrt(sum(im1.^2,1))));