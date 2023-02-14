function JackknifeData = makeJackknifeWallabyOPM(DImageData)
    
    TrialNumber = getTrialNumber(DImageData);
    JackknifeData = zeros([getMapSize(DImageData) TrialNumber]);
    %rng(seed)

    %BootstrapCombinations = randi([1 getTrialNumber(DImageData)],getTrialNumber(DImageData),getStimulusNumber(DImageData),NBootstrapSamples);
%     disp(getMapSize(DImageData))
%     disp(size(BootstrapCombinations))
%     disp([getMapSize(DImageData) 1])
    
    %JackknifeData(:,:,1) = makeWallabyOPMJason(DImageData);

    for i_TrialNumber = 1:TrialNumber
        disp(i_TrialNumber)
        
        disp('make dimg Jackknifesample')
        dimg = getJackknifeSample(DImageData,i_TrialNumber);
        disp('make OPM')
        OPM = makeWallabyOPMJason(dimg);
        disp('store OPM')
        JackknifeData(:,:,i_TrialNumber) = reshape(OPM,[getMapSize(DImageData) 1]);
        
        %AllMaps(:,:,1+i_BootstrapSample) = reshape(makeWallabyOPMJason(getBootstrapSample(DImageData,BootstrapCombinations(:,:,NBootstrapSamples))),[getMapSize(DImageData) 1]);
        
    end
    
    %OPMs = AllMaps(:,:,1);
    
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

function JackknifeSample = getJackknifeSample(ImageData,i_TrialNumber)
    ImageDataArray = cell2array(ImageData);
    JackknifeSet = getJackknifeSet(getTrialNumber(ImageData),i_TrialNumber);
     
%     disp(getStimulusNumber(ImageData))
%     disp(size(JackknifeSet ))
    
    JackknifeSampleArray = zeros([getStimulusNumber(ImageData),getTrialNumber(ImageData)-1,size(ImageDataArray,3:5)]);
    
%     disp([getStimulusNumber(ImageData),size(JackknifeSet ),size(ImageDataArray,3:5)])
%     disp(size(JackknifeSampleArray))
%     
%     disp(JackknifeSet)
    
    %disp(size(ImageDataArray))
    
    for i_Stimulus = 1:getStimulusNumber(ImageData)
        
        %disp(size(BootstrapCombinations(:,i_Stimulus)))
%         disp(size(ImageDataArray(i_Stimulus,JackknifeSet,:,:,:)))
%         disp(size(JackknifeSampleArray(i_Stimulus,:,:,:,:)))
        
        JackknifeSampleArray(i_Stimulus,:,:,:,:) = ImageDataArray(i_Stimulus,JackknifeSet,:,:,:);
    end
    
    JackknifeSample = array2cell(JackknifeSampleArray);
end

function JackknifeSet = getJackknifeSet(TrialNumber,i_TrialNumber)
    JackknifeSet = setdiff(1:TrialNumber,i_TrialNumber);
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
    Delta = abs(AllMaps - AllMaps(:,:,1));
    DeltaAngle = angle(AllMaps) - angle(AllMaps(:,:,1));
    %DeltaMean = mean(abs(Delta),1:2);
    MaxDeltaOPM = AllMaps(:,:,getIndexMaxDeltaOPM(Delta));
    ConfidenceIntervalls(:,:,1) = getConfidenceIntervalls(Delta,Confidence); 
    ConfidenceIntervalls(:,:,2) = acos(-getConfidenceIntervalls(-cos(DeltaAngle),Confidence));
end

function I = getIndexMaxDeltaOPM(Delta)
    [~,I] = max(mean(Delta,1:2));
end

function ConfidenceIntervalls = getConfidenceIntervalls(Delta,Confidence)
    DeltaSort = sort(Delta,3);
    ConfidenceIntervalls = DeltaSort(:,:,round(size(Delta,3)*Confidence));   
end