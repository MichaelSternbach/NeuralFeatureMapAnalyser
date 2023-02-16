function map = make_map(data,stim_order,ROI,do_gif)
% Calculate map from raw data
%
% map = make_map(data,stim_order,ROI,do_gif)
%
% Input:
%  - data is ordered as: [y pixel, x pixel, stimulus condition, block]
%  - stim_order: stimulus conditions. If complex = binocular. NaN = blank.
%  - ROI:  determines the pixels in the image that are part of the map. If
%  left empty, the whole image will be used for averaging.
% Output:
% - map is the complex valued orientation layout, where the magnitude is
% the orientation selectivity and the phase the orientation preference.
%
% chepe@nld.ds.mpg.de
%%

% check input
if ndims(data)<3 || ndims(data)>4
    error('Check data dimensions')
end

[nPixY,nPixX,~,nblocks] = size(data);
if nargin==1
    stim_order = [NaN 0 45 90 135];
    ROI = true(nPixY,nPixX);
    do_gif = false;
elseif nargin==2
    ROI = true(nPixY,nPixX);
    do_gif = false;
elseif nargin==3
    do_gif = false;
end
if isempty(ROI)
    ROI = true(nPixY,nPixX);
end

% apply GIF method to data
if do_gif && nblocks>1
    [data,stim_order] = reduce_data(data,stim_order);
    disp('GIF')
end
    
% average the blocks for each condition
data = mean(data,4);

% make maps for each frame
map = zeros(nPixY,nPixX);
for stim_ii = 1:length(stim_order)
    % skip blanks
    if isnan(real(stim_order(stim_ii)))
        continue
    end
    map = map + data(:,:,stim_ii) * exp(1i*2*pi/180*real(stim_order(stim_ii)));
end

% Z-score the map (common practice)
map = (map-mean(map(ROI)))/std(map(ROI));

end


function [data_clean,stim_unique] = reduce_data(data,stim_order)

[nPixY,nPixX,~,nblocks] = size(data);

% compress data if binocular and direction
stim_unique = unique(mod(real(stim_order),180));
stim_unique(isnan(stim_unique)) = [];

repeats = sum(bsxfun(@minus,stim_unique,mod(real(stim_order),180)')==0,1);
if ~all(repeats==repeats(1))
    % if each stimulus is not repeated in each eye
    data_clean = zeros(nPixY,nPixX,length(stim_unique),nblocks);
    for stim_ii=1:length(stim_unique)
        data_clean(:,:,stim_ii,:) = mean(data(:,:,real(stim_order)==stim_unique(stim_ii),:),3);
    end
else
    data_clean = zeros(nPixY,nPixX,length(stim_unique),nblocks*repeats(1));
    for stim_ii=1:length(stim_unique)
        data_clean(:,:,stim_ii,:) = reshape(...
            data(:,:,mod(real(stim_order),180)==stim_unique(stim_ii),:),...
            [nPixY,nPixX,nblocks*repeats(1)]);
    end
end

data_clean = GIF(data_clean);

end
