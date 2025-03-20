function dataOut = NoAveragePreProcessRawDataJason(expIds,refWin,sigWin,partId,dataPath,ID)
    %addpath '/home/michael/Cloud/PhD/data/data share/Wallaby data/Script/'
    % optical imaging analysis/report

    % 16/9/2014 - Shaun L. Cloherty <s.cloherty@ieee.org>
    % 20/6/2022 - Young J. Jung

    % close all;
    %clear

    % Imaging data is binned 2x2 'in camera' - 1px = 24um.
    %px = 1/0.024; % 1px = 24um

    % reference/baseline window
    %refWin = [1:10]; % average over pre-stimulus images (stimulus onset at 2s)

    % signal window
    %sigWin = [31:35]; % FIXME: average around peak of mean intensity plot?

    %partId = ['A','B']; % 8 orientations
    % expIds = {[3 5 7],[4 6 8]}; %600 microns down
    % expIds = {[15 17 19],[16 18 20]}; %200 microns down
    %expIds = {[29 33 35],[30 34 36]}; %400 microns down
    % % expIds = {[46],[47]}; %Right hemisphere
    % expIds = {[37],[38]};

    smallMem=0;

    %dataPath = '/media/michael/Samsung_T5/Jason/Dunnart Data/DunnartAN';%'D:\Jason\Dunnart Data\DunnartAN';
    expPrefix = '.';

    %grnImg = oiLoadGrnImage( fullfile(dataPath, expPrefix, 'surface3_roi1x1.bmp'));

    % addpath (genpath ('C:\1- PhD work\Data\WallabyH\Optical Imaging;'))
    % dataPath = 'T:\IN VIVO LAB\Data\Dunnart\DunnartAH\Optical Imaging\';% load 'exp%i'
    % dataPath ='T:\IN VIVO LAB\Data\Dunnart\DunnartAH\Optical Imaging';
    % grnImg = oiLoadGrnImage( fullfile(dataPath, 'GrnImg_1620.bmp'));

    %%load images 
    for i = 1:length(partId),
        fprintf(1, 'Part %c: ', partId(i));

        dimg = {};

        for j = 1:length(expIds{i}),
            fname = sprintf([ID 'exp%i.mat'], expIds{i}(j)); %fix this %dunnartAN
            load(fullfile(dataPath, fname));

            fprintf(1, 'exp%i (%ix%i) ', expIds{i}(j), size(images));

    %         yy=[65 220] ;xx =[135 305];% ROI
    %         images = oiCrop(images, yy, xx);

                bimg = images;

            clear images

            sigWin_ = sigWin;
%             if bitand(smallMem,2),
%                 bimg = cellfun(@(x) x(:,:,[refWin,sigWin]), bimg, 'UniformOutput',0);
% 
%             % redefine sigWin (see comment aove)
%             sigWin_ = setdiff(1:size(bimg{1},3),refWin);
%             end

            % Calcuate difference image sequence
            dimg = [dimg, oiDiff(bimg,refWin)];

        end

    %     Combine opposite directions - now looking at 'orientation'
%         [aimg(i:2:8,1)] = oiAve([dimg(1:4,:), dimg(5:8,:)]);
%          aimg_(i,1) = oiAve(dimg(9,:));
        %% Combine opposite directions for each trial
        for i_trial = 1: size(dimg,2)
            [aimg(i:2:8,i_trial)] = oiAve([dimg(1:4,i_trial), dimg(5:8,i_trial)]);
            aimg_(i,i_trial) = oiAve(dimg(9,i_trial));
        end


    % %     Looking at direction
    %     [aimg(1:8,1)] = oiAve(dimg(1:8,:));
    %     aimg_(i,1) = oiAve(dimg(9,:));

    end
    
    
    %% Transfer cells to array
    disp('cell to array')
    data = zeros([size(aimg{1,1},1:2) size(aimg,1)+1 size(aimg,2) size(aimg{1,1},3)]);
    for i_trial = 1: size(aimg,2)
        for i_stim = 1: size(aimg,1)
            data(:,:,i_stim,i_trial,:)=aimg{i_stim,i_trial};
        end
        data(:,:,9,i_trial,:) = (aimg_{1,i_trial}+aimg_{2,i_trial})/2;
    end
    
%     %% Cocktail party applied to aimg
%     disp('Apply coctail party')
%     %data(:,:,1:8,:,:) = data(:,:,1:8,:,:) - mean(data(:,:,1:8,:,:),3);
%     for i = 1: size(data,4)
%         data(:,:,1:8,i,:) = data(:,:,1:8,i,:) - mean(data(:,:,1:8,i,:),3);
%     end
    
    disp('time average')
    dataOut = -mean(data(:,:,:,:,sigWin),5);

    %% remove average
    disp('remove meane')
    dataOut(:,:,1:end-1,:) = dataOut(:,:,1:end-1,:)-mean(dataOut(:,:,1:end-1,:),3:4);

    
end
