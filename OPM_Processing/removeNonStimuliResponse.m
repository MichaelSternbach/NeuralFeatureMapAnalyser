function dataOut = removeNonStimuliResponse(data,NonStimuliResponseRemovalMethod,sigWin)
    
    %% check input andset default values
    if length(size(data)) < 4 || length(size(data)) > 5
        error('data must be 4D or 5D')
    end
    if nargin < 2
        NonStimuliResponseRemovalMethod = 'cocktail party';
    end
    if nargin < 3
        if length(size(data)) == 4
            sigWin = [];
        else
            raise error('sigWin must be provided for 5D data');
        end
    end

    %% time average
    if length(size(data)) == 5 && ~isempty(sigWin)
        disp('time average')
        dataOut = mean(data(:,:,:,:,sigWin),5);
    else
        disp('no time average')
        dataOut = data;
    end


    %% remove non-stimuli response
    switch lower(NonStimuliResponseRemovalMethod)
        case {'average','mean','average response','mean response'}
            disp('remove average')
            dataOut(:,:,1:end-1,:) = dataOut(:,:,1:end-1,:)-nanmean(dataOut,3:4);
        case {'cocktail party','cocktailparty'}
            disp('remove cocktail party')
            for i = 1: size(dataOut,4)
                dataOut(:,:,1:end-1,i) = dataOut(:,:,1:end-1,i) - nanmean(dataOut(:,:,1:end-1,i),3);
            end
        case {'first frame','firstframe','ff'}
            disp('remove first frame')
            if length(size(data)) ==5 
                dataOut(:,:,1:end-1,:) = dataOut(:,:,1:end-1,:)-data(:,:,1:end-1,:,1);
            else
                raise error('first frame removal can only be used with 5D data that includes time')
            end
        case {'none','no','off'}
            disp('do not remove non-stimuli response')
        otherwise
            error('Unknown NonStimuliResponseRemovalMethod method')
    end
end