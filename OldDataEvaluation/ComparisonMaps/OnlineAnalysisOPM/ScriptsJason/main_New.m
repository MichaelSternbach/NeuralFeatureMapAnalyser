% optical imaging analysis/report

% 16/9/2014 - Shaun L. Cloherty <s.cloherty@ieee.org>
% 20/6/2022 - Young J. Jung

close all;
clear

% Imaging data is binned 2x2 'in camera' - 1px = 24um.
px = 1/0.024; % 1px = 24um

% reference/baseline window
refWin = [1:10]; % average over pre-stimulus images (stimulus onset at 2s)

% signal window
sigWin = [31:35]; % FIXME: average around peak of mean intensity plot?

partId = ['A','B']; % 8 orientations
expIds = {[29,31,33,36,38,40,42,44],[30,32,35,37,39,41,43,45]}; 
smallMem=1;

% addpath (genpath ('C:\1- PhD work\Data\WallabyH\Optical Imaging;'))
dataPath = 'C:\NewMatlabdata\';% load 'exp%i'

%%load images 
for i = 1:length(partId),
    fprintf(1, 'Part %c: ', partId(i));
    
    dimg = {};
    
    for j = 1:length(expIds{i}),
        fname = sprintf('wallabyHexp%i.mat', expIds{i}(j));
        load(fullfile(dataPath, fname));
        
        fprintf(1, 'exp%i (%ix%i) ', expIds{i}(j), size(images));
        
        yy=[65 220] ;xx =[135 305];% ROI
        images = oiCrop(images, yy, xx);

            bimg = images;

        clear images
       
        sigWin_ = sigWin;
        if bitand(smallMem,2),
            bimg = cellfun(@(x) x(:,:,[refWin,sigWin]), bimg, 'UniformOutput',0);
            
        % redefine sigWin (see comment aove)
        sigWin_ = setdiff(1:size(bimg{1},3),refWin);
        end
        
        % Calcuate difference image sequence
        dimg = [dimg, oiDiff(bimg,refWin)];
        
    end
    
%     Combine opposite directions - now looking at 'orientation'
    [aimg(i:2:8,1)] = oiAve([dimg(1:4,:), dimg(5:8,:)]);
     aimg_(i,1) = oiAve(dimg(9,:));
     
%%     Looking at direction
%     [aimg(i:2:16,1)] = oiAve(dimg(1:8,:));
%     aimg_(i,1) = oiAve(dimg(9,:));
%     
end

%% Cocktail party applied to aimg
aimg_sum =(aimg{1,1}+aimg{2,1}+aimg{3,1}+aimg{4,1}+aimg{5,1}+aimg{6,1}+aimg{7,1}+aimg{8,1})./8;
    for j = 1:size(aimg,1),
        caimg{j,1}= aimg{j,1}-aimg_sum;
        
    end

%% High pass filter average image sequence
for i = 1:8
    for k =1:50
        sigma = 20;
        h = fspecial('gaussian', 6*sigma, sigma);
        temp1 = imfilter(caimg{i,1}(:,:,k),h,'symmetric','conv');  
        faimg{i,1}(:,:,k) =caimg{i,1}(:,:,k)-temp1; 
    end
end

%% Low pass filter average image sequence
for i = 1:8
    for k =1:50
        sigma = 2; 
        h = fspecial('gaussian', 6*sigma, sigma);
       faimg{i,1}(:,:,k) = imfilter(faimg{i,1}(:,:,k),h,'symmetric','conv');  
       
    end
end

%% Extended Spatial Decorrelation Method 
dr = [5, 5];
c=[];

for i= 1:size(faimg,1)
    [simg{i,1}(:,:,:), c{i,1}] = esd(faimg{i,1}, dr);
end

%% Calculate single condition maps
numConds = size(faimg,1);
for j = 1:numConds,
    maps(:,:,j) = -1*mean(faimg{j,1}(:,:,sigWin_),3);
end

%% Plot coefficents over time series 
figure;
for j=1:8
    
for i=[1:5,45:50];
    W=real(inv(c{j,1}));
    bleft = maps(:,:,j);
    bright= real(simg{j,1}(:,:,i));
    F=W(:,i)*sign(bleft(:)'*bright(:));
   
    subplot(2,4,j)
    plot(F)
    hold on;
    
end
end

%% coefficient
W =[];
for i=1:8
 W{i,1}=inv(c{i,1});
end

%% Generate maps based on ESD % Fix me
maps(:,:,1) =(-simg{1,1}(:,:,1))*mean(W{1,1}(31:35,1),1);
maps(:,:,2) =(-simg{2,1}(:,:,1))*mean(W{2,1}(31:35,1),1);
maps(:,:,3) =(-simg{3,1}(:,:,1))*mean(W{3,1}(31:35,1),1);
maps(:,:,4) =(-simg{4,1}(:,:,50))*mean(W{4,1}(31:35,50),1);
maps(:,:,5) =(-simg{5,1}(:,:,1))*mean(W{5,1}(31:35,1),1);
maps(:,:,6) =(-simg{6,1}(:,:,50))*mean(W{6,1}(31:35,50),1);
maps(:,:,7) =(-simg{7,1}(:,:,1))*mean(W{7,1}(31:35,1),1);
maps(:,:,8) =(-simg{8,1}(:,:,50))*mean(W{8,1}(31:35,50),1);

%% low-pass filter to remove high-frequency noise
sigma = 8; %Fix me 
h = fspecial('gaussian', 6*sigma, sigma);
for i=1:8
fmaps(:,:,i) = imfilter(maps(:,:,i),h,'symmetric','conv');
end

%% Plot OP color maps

figure(1);
[op,r] = oiCalcORMap(real(fmaps), [0:7]*pi/8, 'vector'); % 8 orientations
op = kron(op,ones([1,1])); r = kron(r,ones([1,1]));
h= subimage(uint8(256*op/pi), hsv(256));
set(h,'AlphaData', min(sqrt(r/max(r(:)/6)), 1));
set(gca,'xtick',[]); set(gca,'ytick',[]);


%% Calculate pinwheels on orientation maps
angles = [-90:11.25:78.75];
op=op-degtorad(90);
temp = double(op);
map_pref = temp*180/pi;

[pinw, windno] = locate_pinwheels(map_pref);

% Calculate pinwheel density and column spacing
[~, density] = pinw_density(map_pref);

%% Plot orientation preference distribution 

figure(2);
oiDistribution(real(fmaps));

