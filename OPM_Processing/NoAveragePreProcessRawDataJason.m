function dataOut = NoAveragePreProcessRawDataJason(expIds,refWin,sigWin,partId,dataPath,ID,getDirectionData,NonStimuliResponseRemoval)
    
    if nargin < 7
        getDirectionData = false;
    end

    if nargin <8
        % NonStimuliResponseRemoval = 'average';
        NonStimuliResponseRemoval = 'cocktail party';
    end



 
    for i = 1:length(partId)
        fprintf(1, 'Part %c: ', partId(i));

        dimg = {};

        for j = 1:length(expIds{i})
            fname = sprintf([ID 'exp%i.mat'], expIds{i}(j)); %fix this %dunnartAN
            load(fullfile(dataPath, fname));

            fprintf(1, 'exp%i (%ix%i) ', expIds{i}(j), size(images));

            bimg = images;

            clear images

            % Calcuate difference image sequence
            dimg = [dimg, oiDiff(bimg,refWin)];

        end

    %     Combine opposite directions - now looking at 'orientation'
%         [aimg(i:2:8,1)] = oiAve([dimg(1:4,:), dimg(5:8,:)]);
%          aimg_(i,1) = oiAve(dimg(9,:));
        
        %% loop over trials
        for i_trial = 1: size(dimg,2)
        
            if ~getDirectionData
                % Combine opposite directions - now looking at 'orientation'
                [aimg(i:2:8,i_trial)] = oiAve([dimg(1:4,i_trial), dimg(5:8,i_trial)]);
                aimg_(i,i_trial) = oiAve(dimg(9,i_trial));
            else
                % Looking at direction
                [aimg(i:2:16,i_trial)] = oiAve(dimg(1:8,i_trial));
                aimg_(i,i_trial) = oiAve(dimg(9,i_trial));
            end
        end




    end
    
    
    %% Transfer cells to array
    disp('cell to array')
    data = zeros([size(aimg{1,1},1:2) size(aimg,1)+1 size(aimg,2) size(aimg{1,1},3)]);
    for i_trial = 1: size(aimg,2)
        for i_stim = 1: size(aimg,1)
            data(:,:,i_stim,i_trial,:)=aimg{i_stim,i_trial};
        end
        data(:,:,end,i_trial,:) = (aimg_{1,i_trial}+aimg_{2,i_trial})/2;
    end
    
%     %% Cocktail party applied to aimg
%     disp('Apply coctail party')
%     %data(:,:,1:8,:,:) = data(:,:,1:8,:,:) - mean(data(:,:,1:8,:,:),3);
%     for i = 1: size(data,4)
%         data(:,:,1:8,i,:) = data(:,:,1:8,i,:) - mean(data(:,:,1:8,i,:),3);
%     end

    dataOut = processTimeSeriesData(data,sigWin,NonStimuliResponseRemoval);
    

    % %% remove average
    % disp('remove average')
    % dataOut(:,:,1:end-1,:) = dataOut(:,:,1:end-1,:)-mean(dataOut(:,:,1:end-1,:),3:4);

    
end

function dataOut = processTimeSeriesData(data,sigWin,NonStimuliResponseRemoval)
    
    %% time average
    disp('time average')
    dataOut = -mean(data(:,:,:,:,sigWin),5);


    %% remove non-stimuli response
    switch lower(NonStimuliResponseRemoval)
        case {'average','mean','average response','mean response'}
            disp('remove average')
            dataOut(:,:,1:end-1,:) = dataOut(:,:,1:end-1,:)-mean(dataOut,3:4);
        case {'cocktail party','cocktailparty'}
            disp('remove cocktail party')
            for i = 1: size(dataOut,4)
                dataOut(:,:,1:end-1,i) = dataOut(:,:,1:end-1,i) - mean(dataOut(:,:,1:end-1,i),3);
            end
        case {'first frame','firstframe','ff'}
            disp('remove first frame')
            dataOut(:,:,1:end-1,:) = dataOut(:,:,1:end-1,:)+data(:,:,1:end-1,:,1);
        otherwise
            error('Unknown NonStimuliResponseRemoval method')
    end
end