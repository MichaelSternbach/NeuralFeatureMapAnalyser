function maps = makeWallabyOPM(ImageData,StimulusOnset)
addpath '/home/michael/Cloud/git/vone/MatlabCode/PatchAnalysis'
addpath '/home/michael/Cloud/PhD/data/data share/Wallaby data/Script'

%'/home/michael/Cloud/PhD/data/data share/Wallaby data/Maps/WallabyH_Image.mat'
dr = [5, 5];
%plot for testing
StimPlot=3;
TimePlot = 40;


%img = subtractStimulusIndependentPatterns(loadImageData(ImageDataPath),StimulusOnset);

img = subtractStimulusIndependentPatterns(ImageData,StimulusOnset);

faimg = averageTrials(img);
plotGrayMap(reshape(faimg(StimPlot,:,:,TimePlot),size(faimg,2),size(faimg,3)))

faimg = faimg(1:8,:,:,:);
plotGrayMap(reshape(mean(faimg(1,:,:,TimePlot),1),size(faimg,2),size(faimg,3)))

faimg = removeCockTailBlank(faimg);
plotGrayMap(reshape(faimg(StimPlot,:,:,TimePlot),size(faimg,2),size(faimg,3)))


faimg = HighPassFilter(faimg);
plotGrayMap(reshape(faimg(StimPlot,:,:,TimePlot),size(faimg,2),size(faimg,3)))


faimg = LowPassFilter(faimg);
plotGrayMap(reshape(faimg(StimPlot,:,:,TimePlot),size(faimg,2),size(faimg,3)))


SourceFrames = applyESD(faimg,StimulusOnset,dr);
plotGrayMap(reshape(SourceFrames(1,:,:,TimePlot),size(SourceFrames,2),size(SourceFrames,3)))


maps = AveragesPeakIntensityFrames(SourceFrames,5);
plotGrayMap(reshape(maps(:,:,StimPlot),size(maps,1),size(maps,2)))


maps = LowPassFilterMaps(maps);
plotGrayMap(reshape(maps(:,:,StimPlot),size(maps,1),size(maps,2)))

OPM = makeOPM(maps(1:8,:,:));

figure()
plot_map(OPM)
disp(size(OPM))
plotOPM(maps)


end


function dimg = loadImageData(ImageDataPath)
    load(ImageDataPath);
end

function img = subtractStimulusIndependentPatterns(img,StimulusOnset)
    %disp(size(img)) 
    for i_stimulus = 1:size(img,1)
       for i_trial = 1:size(img,2)
          img(i_stimulus,i_trial) =  {subtractStimulusIndependentPattern(img(i_stimulus,i_trial),StimulusOnset)};
       end
    end
end

function img_cell = subtractStimulusIndependentPattern(cell,StimulusOnset)
    img_cell = cell{1}-mean(cell{1}(:,:,1:StimulusOnset),3);
end

function faimg = averageTrials(img)
    DimFrames = size(oc(img));
    faimg = zeros(size(img,1),DimFrames(1),DimFrames(2),DimFrames(3),'single');
    for i_stimulus = 1:size(img,1)
        for i_trial = 1:size(img,2)
          faimg(i_stimulus,:,:,:)=faimg(i_stimulus,:,:,:)+reshape(oc(img(i_stimulus,i_trial)),1,DimFrames(1),DimFrames(2),DimFrames(3));
       end
    end
    faimg = faimg./size(img,2);
end

function opend_cell = oc(cell)
    opend_cell = cell{1};
end

function faimg = removeCockTailBlank(faimg)
    faimg = faimg - mean(faimg,1);
end

function faimg = HighPassFilter(faimg)
    sigma = 20;
    h = fspecial('gaussian', 6*sigma, sigma);
    for i_stimulus = 1:size(faimg,1)
        for i_frame = 1:size(faimg,4)
          faimg(i_stimulus,:,:,i_frame)=faimg(i_stimulus,:,:,i_frame)-imfilter(faimg(i_stimulus,:,:,i_frame),h,'symmetric','conv'); 
       end
    end
end

function faimg = LowPassFilter(faimg)
    sigma = 2;
    h = fspecial('gaussian', 6*sigma, sigma);
    for i_stimulus = 1:size(faimg,1)
        for i_frame = 1:size(faimg,4)
          faimg(i_stimulus,:,:,i_frame)=imfilter(faimg(i_stimulus,:,:,i_frame),h,'symmetric','conv'); 
       end
    end
end

function SourceFrames = applyESD(faimg,StimulusOnset,dr)
    
    SourceFrames = zeros(size(faimg));
    for i_stimulus = 1:size(faimg,1)
        %disp(size(faimg(i_stimulus,:,:,:)))
        [Y, W] = esd(reshape(faimg(i_stimulus,:,:,:),size(faimg,2),size(faimg,3),size(faimg,4)),dr);
        SourceFrames(i_stimulus,:,:,:) = makeSourcesFromESD(Y, W,StimulusOnset);
    end
end

function SourceFrames = makeSourcesFromESD(Y, W,StimulusOnset)
%     SourceFrames = zeros(size(Y))
%     for i_source = 1:size(Y,3)
%         SourceFrames(:,:,i_sources) =simg(:,:,i_sources)*inv(W(:,i_sources);
%     end
%[Y, W] = esd(reshape(faimg(i_stimulus,:,:,:),size(faimg,2),size(faimg,3),size(faimg,4)),dr);
    B=inv(W);
    StimulusSource = maxSignalStimulus(B,StimulusOnset);
    PlotSourceStrengthes(B,StimulusSource);
    SourceFrames = reshape(Y(:,:,StimulusSource),size(Y,1),size(Y,2),1).*reshape(B(:,StimulusSource),1,1,size(B,1));
end


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


function PlotSourceStrengthes(W,StimulusSource)
    figure()
    %title(TitleText)
    %title('test')
    for i_source = [StimulusSource 1:6]
       %disp(size(W(:,i_source)))
       plot(1:size(W,1),real(reshape(W(:,i_source),1,size(W,1))),'DisplayName',num2str(i_source))%reshape(W(:,i_source),size(W,1))
       hold on
    end
    hold off
    legend
end

function StimulusSource = maxSignalStimulus(W,StimulusOnset)
    [~,StimulusSource]=max(getSignalStimulus(W,StimulusOnset));
end

function SignalStimulus= getSignalStimulus(W,StimulusOnset)
    SignalStimulus = std(W(StimulusOnset:size(W,1),:),1);
end

% function maps = AveragesPeakIntensityFrames(SourceFrames,NAverage) 
%     maps = zeros(size(SourceFrames,1),size(SourceFrames,2),size(SourceFrames,3),'single');
%     for i_stimulus = 1:size(SourceFrames,1)
%         maps(i_stimulus,:,:) = averagePeakIntensityFrames(SourceFrames(i_stimulus,:,:,:),NAverage) ;
% 
%     end
% end


function maps = AveragesPeakIntensityFrames(SourceFrames,NAverage) 
    maps = zeros(size(SourceFrames,2),size(SourceFrames,3),size(SourceFrames,1),'single');
    for i_stimulus = 1:size(SourceFrames,1)
        maps(:,:,i_stimulus) = reshape(averagePeakIntensityFrames(SourceFrames(i_stimulus,:,:,:),NAverage),size(SourceFrames,2),size(SourceFrames,3),1);

    end
end

function m = averagePeakIntensityFrames(SourceFrames,NAverage) 
    %peak = mean(getPeak(SourceFrames));
    lb=35;%31; %(peak-(NAverage-1)/2);
    ub=40;%35; %(peak+(NAverage-1)/2);
    %disp(size(SourceFrames))
    %disp(peak)
    m = mean(SourceFrames(:,:,:,lb:ub),4);
end

function peak = getPeak(SourceFrames)
    [~,peak] = max(mean(mean(abs(SourceFrames),2),3));
end

function maps = LowPassFilterMaps(maps)
    sigma = 8; 
    h = fspecial('gaussian', 6*sigma, sigma);
    for i_stimulus = 1:size(maps,3)
        maps(:,:,i_stimulus) = imfilter(maps(:,:,i_stimulus),h,'symmetric','conv');

    end
end 

% function OPM = makeOPM(maps)
% OPM = reshape(mean(maps.*exp(2i*pi/size(maps,3)*reshape([1:size(maps,3)],1,1,size(maps,3))),3),size(maps,1),size(maps,2));
% end
function plotGrayMap(map)
    figure()
    a = abs(map);
    %figure;
    imagesc(a); 
    %colormap jet;
    %colormap hsv;
    colormap gray;
end

function plotOPM(fmaps)
    figure()
    [op,r] = oiCalcORMap(real(fmaps), [0:7]*pi/8, 'vector'); % 8 orientations
    op = kron(op,ones([1,1])); r = kron(r,ones([1,1]));
    h= imshow(uint8(256*op/pi), hsv(256));
    set(h,'AlphaData', min(sqrt(r/max(r(:)/6)), 1));
    set(gca,'xtick',[]); set(gca,'ytick',[]);
end
