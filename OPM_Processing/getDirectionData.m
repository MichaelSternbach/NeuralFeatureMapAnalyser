function [DirectionMap,StimDir] = getDirectionData(data_info,data_path,GIF_SN_TH)
    
    SpacingDirectionData = [data_path '/DirectionData_' data_info.animal data_info.ID '.mat'];
    if isfile(SpacingDirectionData)
        load(SpacingDirectionData,'DirectionData')
    else
        DirectionData = NoAveragePreProcessRawDirectionDataJason(data_info.expIds,data_info.refWin,data_info.sigWin,data_info.partId,data_path,data_info.ID);
        save(SpacingDirectionData,'DirectionData')
    end
    
    if GIF_SN_TH >0
        DirectionData = GIF(DirectionData,GIF_SN_TH);
    end

    StimDir = [data_info.stim_order(find(~isnan(data_info.stim_order))) data_info.stim_order(find(~isnan(data_info.stim_order)))+180];
    StimDir = reshape(StimDir,[1 1 length(StimDir) 1]);
    DirectionMap = mean(DirectionData(:,:,1:size(DirectionData,3)-1,:).*exp(2i*pi*StimDir/360),3:4);
end