function [AllMaps,MaXDeltaOPM,ConfidenceIntervalls] = VariabilityTestWallabyOPM(DImageData,NBootstrapSamples,seed)
    AllMaps = zeros([getMapSize(DImageData) NBootstrapSamples+1]);
    rng(seed)

    BootstrapCombinations = randi([1 getTrialNumber(DImageData)],getTrialNumber(DImageData),getStimulusNumber(DImageData),NBootstrapSamples);
    
%     disp(getMapSize(DImageData))
%     disp(size(BootstrapCombinations))
%     disp([getMapSize(DImageData) 1])
    
    AllMaps(:,:,1) = makeWallabyOPMJason(DImageData);

    for i_BootstrapSample = 1:NBootstrapSamples
        disp(i_BootstrapSample)
        
        disp('make dimg Bootstrapsample')
        dimg = getBootstrapSample(DImageData,BootstrapCombinations(:,:,i_BootstrapSample));
        disp('make OPM')
        OPM = makeWallabyOPMJason(dimg);
        disp('store OPM')
        AllMaps(:,:,1+i_BootstrapSample) = reshape(OPM,[getMapSize(DImageData) 1]);
        
        %AllMaps(:,:,1+i_BootstrapSample) = reshape(makeWallabyOPMJason(getBootstrapSample(DImageData,BootstrapCombinations(:,:,NBootstrapSamples))),[getMapSize(DImageData) 1]);
        
    end
    
    %OPMs = AllMaps(:,:,1);
    [MaXDeltaOPM,ConfidenceIntervalls]= getVariability(AllMaps);
    
end

function MapSize = getMapSize(ImageData)
    MapSize = size(ImageData{1,1},1:2);
end

function TrialNumber = getTrialNumber(ImageData)
    TrialNumber = size(ImageData,2);
end

function FrameNumber = getFrameNumber(ImageData)
    FrameNumber = size(ImageData{1,1},3);
end

function StimulusNumber = getStimulusNumber(ImageData)
    StimulusNumber = size(ImageData,1);
end

function BootstrapSample = getBootstrapSample(ImageData,BootstrapCombinations)
    ImageDataArray = cell2array(ImageData);
    BootstrapSampleArray = ImageDataArray;
          
    %disp(size(ImageDataArray))
    for i_Stimulus = 1:getStimulusNumber(ImageData)
        %disp(size(BootstrapCombinations(:,i_Stimulus)))
        BootstrapSampleArray(i_Stimulus,:,:,:,:) = ImageDataArray(i_Stimulus,BootstrapCombinations(:,i_Stimulus),:,:,:);
    end
    
    BootstrapSample = array2cell(BootstrapSampleArray);
end

function array = cell2array(cell)
    array = zeros([size(cell) size(cell{1,1})]);
    for i_x = 1: size(cell,1)
       for i_y = 1:size(cell,2)
           array(i_x,i_y,:,:,:)=reshape(cell{i_x,i_y},[1 1 size(cell{i_x,i_y})]);
       end
    end
end


function Cell = array2cell(array)
    Cell=cell(size(array,1),size(array,2));
    for i_x = 1: size(array,1)
       for i_y = 1:size(array,2)
           Cell{i_x,i_y} = reshape(array(i_x,i_y,:,:),size(array,3:5));
       end
    end
end

function [MaxDeltaOPM,ConfidenceIntervalls]= getVariability(AllMaps)
    Confidence = 0.9;
    Delta = AllMaps - AllMaps(:,:,1);
    %DeltaMean = mean(abs(Delta),1:2);
    MaxDeltaOPM = max(abs(Delta),[],3);
    ConfidenceIntervalls = getConfidenceIntervalls(Delta,Confidence);
end

function ConfidenceIntervalls = getConfidenceIntervalls(Delta,Confidence)
    DeltaSort = sort(abs(Delta),3);
    ConfidenceIntervalls = DeltaSort(:,:,round(size(Delta,3)*Confidence));
end