function [data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,experiment_num)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [data_info,data_path] = info_handle(animal,experiment_num);
    data_info.pixels_per_mm = data_info.pix_per_mm;

    switch lower(animal)
        case{'dunnart'}
            data = NoAveragePreProcessRawDataJason(data_info.expIds,data_info.refWin,data_info.sigWin,data_info.partId,data_path,data_info.ID);
            data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
            BloodVesselImg = getBloodVesselImgMarsupial(data_path,data_obj.info.ID);
        case{'ferret'}
            set_blocks = data_info.protocol.blocks;
            trials_to_use = find(set_blocks>0);
            if sum(set_blocks>0,'all')>1
                disp('Trials to use:')
                disp(trials_to_use)
                trial = input('Input trial to use for plots!');
            end
            load([data_path,'/Processed_2/trial_',num2str(trial),'.mat'],'data')
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order);
            data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
        otherwise
            error('animal entry not recognized')
    end

    

end

function BloodVesselImg = getBloodVesselImgFromNanStim(data,stim_order)
    NanStim = find(isnan(stim_order));
    BloodVesselImg = mean(data(:,:,NanStim,:),4);
end