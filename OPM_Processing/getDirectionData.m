function [DirectionMap,StimDir,DirectionData] = getDirectionData(data_info,data_path,GIF_SN_TH,NonStimuliResponseRemoval)
    if nargin < 4
        NonStimuliResponseRemoval = 'cocktail party';
    end
    
    NonStimuliResponseRemoval_txt = replace(NonStimuliResponseRemoval,' ','_');

    SpacingDirectionData = [data_path 'DirectionData_' data_info.animal data_info.ID '_' NonStimuliResponseRemoval_txt '.mat'];
    disp(isfile(SpacingDirectionData))
    disp(SpacingDirectionData)
    if isfile(SpacingDirectionData)
        load(SpacingDirectionData,'DirectionData')
    else
        getDirectionData = true;
        DirectionData = NoAveragePreProcessRawDataJason(data_info.expIds,data_info.refWin,data_info.sigWin,data_info.partId,data_path,data_info.ID,getDirectionData,NonStimuliResponseRemoval);
        save(SpacingDirectionData,'DirectionData')
    end
    
    if GIF_SN_TH >0
        DirectionData(:,:,1:end-1,:) = GIF(DirectionData(:,:,1:end-1,:),GIF_SN_TH);
    end

    StimDir = [data_info.stim_order(find(~isnan(data_info.stim_order))) data_info.stim_order(find(~isnan(data_info.stim_order)))+180];
%     StimDir = [StimDir NaN];
    StimDir = reshape(StimDir,[1 1 length(StimDir) 1]);
    DirectionMap = mean(DirectionData(:,:,1:size(DirectionData,3)-1,:).*exp(2i*pi*StimDir/360),3:4)+pi/2;
end