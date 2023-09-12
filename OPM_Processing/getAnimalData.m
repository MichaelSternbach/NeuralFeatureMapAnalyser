function [data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,experiment_num,trial_ii,DoRectangleROI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
%     %% connect to data server
%     command = 'sshfs sternbach1@login-dbn02.hpc.gwdg.de:/home/uni08/cidbn1/ ~/CIDBN1 -o IdentityFile=~/ssh_key_gwdg_hpc.pu';
%     [status,cmdout] = system(command);
%     disp([status,cmdout])
    
    %% get data infos
    [data_info,data_path] = info_handle(animal,experiment_num);
    if isfield(data_info,'pix_per_mm')
        data_info.pixels_per_mm = data_info.pix_per_mm;
    else
        data_info.pix_per_mm = data_info.pixels_per_mm;
    end

    switch lower(animal)
        case{'dunnart'}
            %% load data
            ProcessedDataFile = [data_path 'ProcessedData.mat']; 
            if isfile(ProcessedDataFile)
                load(ProcessedDataFile,'data')
            else
                data = NoAveragePreProcessRawDataJason(data_info.expIds,data_info.refWin,data_info.sigWin,data_info.partId,data_path,data_info.ID);
                save(ProcessedDataFile,'data')
            end
            
            %% update pix_per_mm
%             data_info.pix_per_mm = 20/0.48*2;
%             data_info.pixels_per_mm = data_info.pix_per_mm;
%             save_info(data_path,data_info)
            
            %% get ROI
            info_path = [data_path,'exp_info.mat'];
            tmp = load(info_path,'ROI');
                if isfield(tmp,'ROI')
                    ROI = tmp.ROI;
                else
                    tmp2 = load([data_path data_info.ID '_Mask.mat'],'a');
                    if isfield(tmp2,'a')
                        ROI = tmp2.a;
                        save(info_path,'ROI',"-append")
                    else
                        disp('ROI is not given!')
                        ROI = true(size(data,1),size(data,2));
                    end
                end
            
            %% make data object
            data_obj = data_handle_corrected(data_info,data,ROI);
            
            %% load blodvessel image
            BloodVesselImg = getBloodVesselImgMarsupial(data_path,data_obj.info.ID);
                        
        case{'ferret'}
            
            %% get trials to use
            set_blocks = data_info.protocol.blocks;
            trials_to_use = find(set_blocks>0);
            if nargin == 2 || trial_ii == 0
                if sum(set_blocks>0,'all')>1
                    disp('Trials to use:')
                    disp(trials_to_use)
                    trial = input('Input trial to use for plots!');
                end
            else
                trial = trials_to_use(trial_ii);
            end
            
            %% load data
            load([data_path,'/Processed_2/trial_',num2str(trial),'.mat'],'data')
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order);
            
            %% make data object
            data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
        case{'macaque_sam'}
             %% load data
            load([data_path,'Processed/trial_1_corrected.mat'],'data')
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order);
            
            %% make data object
            data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
        case {'microcebus','mouse lemur'}
             if nargin == 2 || trial_ii == 0
                disp('Data to use:')
                disp(size(data_info.exp_data,2))
                trial = input('Input which datafolder to use!');
            else
                trial = trial_ii;
            end
            
            %% load data
            DataFile = [data_path 'Processed/' data_info.exp_data(trial).folder data_info.exp_data(trial).binocular{1}];
            load(DataFile,'data')
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order.binocular);
            
            %% make data object
            data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
        case {'macaque_sam'}
            
            %% load data
            DataFile = [data_path 'Processed/trial_1_corrected.mat'];
            load(DataFile,'data')
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order.binocular);
            
            %% make data object
            data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
            
        otherwise
            error('animal entry not recognized')
    end

    %% get rectangle
    if nargin <=3
        DoRectangleROI = false;
    end
    if (~isfield(data_info, 'rectangle') && DoRectangleROI) || DoRectangleROI == -1
        z_input = data_obj.filter_map(data_obj.read_map(1));
        figure
        plot_map(z_input,0,1,ROI)
        rectangle = input('Input a fitting rectangle to cut out from the OPM! Format [xUpLeft yUpLeft xDownRight yDownRight]');
        close
        data_info.rectangle=rectangle;
        save_info(data_path,data_info)
    end
    
    if isfield(data_info, 'rectangle') 
        data_info.RectangleROI = getRectangleROI(data_info.rectangle,ROI);
    else
        data_info.rectangle=[0 0 size(data,1) size(data,2)];
        data_info.RectangleROI = data_obj.ROI;
    end

end

function BloodVesselImg = getBloodVesselImgFromNanStim(data,stim_order)
    NanStim = find(isnan(stim_order));
    BloodVesselImg = mean(data(:,:,NanStim,:),4);
end

function ROI = getROI(data_path,ID)
    load([data_path ID '_Mask.mat'],'a')
    ROI = a;
end