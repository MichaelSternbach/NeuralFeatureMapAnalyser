function [data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,experiment_num,AnimalDataFolder,data_info_file_list,trial_ii,DoRectangleROI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
%     %% connect to data server
%     command = 'sshfs sternbach1@login-dbn02.hpc.gwdg.de:/home/uni08/cidbn1/ ~/CIDBN1 -o IdentityFile=~/ssh_key_gwdg_hpc.pu';
%     [status,cmdout] = system(command);
%     disp([status,cmdout])
    if nargin < 4 %|| isempty(data_info_file_list)
        function_dir = fileparts(mfilename('fullpath')); % Get function directory
        data_info_file_list = fullfile(function_dir, 'experiment_info.csv');
    end

    if nargin <4
        AnimalDataFolder = '~/CIDBN/';
    end

    if nargin<5
        trial_ii =1;
    end
    
    if nargin<6
        DoRectangleROI = false;
    end
    
    animal = strrep(animal,'-',' ');
    
    %% get data infos
    [data_info,data_path] = info_handle(animal,experiment_num,AnimalDataFolder,data_info_file_list);
    if isfield(data_info,'pix_per_mm')
        data_info.pixels_per_mm = data_info.pix_per_mm;
    else
        data_info.pix_per_mm = data_info.pixels_per_mm;
    end
    
    %% load data in different ways based on the type of animal 
    switch lower(animal)

        case{'mouse'}
            %% load data
            load([data_path,data_info.ID,'.mat'],'data')

            %% change stim order
            data_info.stim_order = [0 45 90 135];

            %% convert Kalatsky data
            data = convertKalatskyData(data,data_info);
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromMean(data);

            %% load ROI
            load([data_path,data_info.ID,'.mat'],'corrections')
%             ROI = corrections.ROI_inside+corrections.ROI_outside;
            ROI = ones(size(data,[1 2]));
            
            %% make data object
            data_obj = data_handle_corrected(data_info,data,ROI);

        case{'dunnart','cat_jung'}
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
                    MaskFile = [data_path data_info.ID '_Mask.mat'];
                    if isfile(MaskFile)
                        tmp2 = load(MaskFile,'a');
                    else
                        tmp2 = struct();
                    end
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
            try 
                BloodVesselImg = getBloodVesselImgMarsupial(data_path,data_obj.info.ID);
            catch
                BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order);
            end
            
            
        case{'wallaby'}
            %% load data
            ProcessedDataFile = [data_path 'ProcessedData.mat']; 
            if isfile(ProcessedDataFile)
                load(ProcessedDataFile,'data')
            else
                data = NoAveragePreProcessRawDataJason(data_info.expIds,data_info.refWin,data_info.sigWin,data_info.partId,data_path,data_info.ID);
                save(ProcessedDataFile,'data')
            end
            
            %% get field_size_pix
            data_info.field_size_pix = size(data,[1,2]);
            
            %% make data object
            ROI = true(size(data,1),size(data,2));
            data_obj = data_handle_corrected(data_info,data,ROI);
            
            %% get new ROI
            %ROI = true(size(data,1),size(data,2));%
            %ROI = getRectangleROI(data_info.rectangle,true(size(data,1),size(data,2)));
            ROI_File = [data_path animal '_SetROI_' data_obj.info.ID];
            load([ROI_File '.mat'],'ROI')
            data_obj.set_ROI(ROI)
            info_path = [data_path,'exp_info.mat'];
            save(info_path,'ROI',"-append")

            
%             info_path = [data_path,'exp_info.mat'];
%             tmp = load(info_path,'ROI');
%             if isfield(tmp,'ROI')
%                 ROI = tmp.ROI;
%             else
%                 disp('calc ROI via column spacing')
%                 smallest_w_mm = 0.01;
%                 w_step_mm = 0.05;
%                 largest_w_mm = 1.5;
%                 [~,~,newROI,~] = get_column_spacingManuel(data_obj.filter_map(data_obj.read_map()),data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
%                 if ~any(newROI)
%                     ROI = newROI;   
%                 else
%                     disp('ROI is not given!')
%                     ROI = true(size(data,1),size(data,2));
%                 end
%                 data_obj.set_ROI(ROI)
%                 save(info_path,'ROI',"-append")
%             end
            
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order);
%             %% load blodvessel image
%             BloodVesselImg = getBloodVesselImgMarsupial(data_path,data_obj.info.ID);
            
                        
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
            %data_obj.apply_LSM()


        case{'cat'}
            
            %% get trials to use
            if nargin == 2 || trial_ii == 0
                TriaToUse =trial_ii;
            else
                TriaToUse =1;
            end
            
            %% load data
            load([data_path,'/PreProcessed/trial_',num2str(TriaToUse),'.mat'],'data')
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order);
            
            %% ROI
            ROI = ones(size(data,1:2));

            %% make data object
            data_obj = data_handle_corrected(data_info,data,ROI);
            %data_obj.apply_LSM()
        case{'cats_all'}
            
            %% get trials to use
            if nargin == 2 || trial_ii == 0
                TriaToUse =trial_ii;
            else
                TriaToUse =1;
            end
            
            %% load data
            load([data_path,'trial_',num2str(TriaToUse),'.mat'],'data')
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order);
            
            %% ROI
            ROI = ones(size(data,1:2));

            %% make data object
            data_obj = data_handle_corrected(data_info,data,ROI);
            %data_obj.apply_LSM()

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
            %data_obj.apply_LSM()
        case {'macaque_sam'}
            
            %% load data
            DataFile = [data_path 'Processed/trial_1_corrected.mat'];
            load(DataFile,'data')
            
            %% make blodvessel image
            BloodVesselImg = getBloodVesselImgFromNanStim(data,data_info.stim_order);
            
            %% make data object
            data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
            %data_obj.apply_LSM()
        otherwise
            error('animal entry not recognized')
    end
    
    %% save data_path
    data_obj.set_data_path(data_path);

    %% set field_size_pix
    if ~isfield(data_info,'field_size_pix')
        data_info.field_size_pix = size(data,1:2);
    end

    %% get rectangle
    if nargin <=3
        DoRectangleROI = false;
    end
    if size(DoRectangleROI,2) == 4
        
        data_info.rectangle=DoRectangleROI;       
        data_info.RectangleROI = getRectangleROI(data_info.rectangle,data_obj.ROI);
        data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
        
    elseif DoRectangleROI == 'max'
        
        [YROI,XROI] = find(data_obj.ROI);
        [Xmin, Xmax] = findBorders(XROI);
        [Ymin, Ymax] = findBorders(YROI);
        
        data_info.rectangle=[Xmin Ymin Xmax Ymax];       
        data_info.RectangleROI = getRectangleROI(data_info.rectangle,data_obj.ROI);
        data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
        
    elseif (~isfield(data_info, 'rectangle') && DoRectangleROI) || DoRectangleROI == -1
        z_input = data_obj.filter_map(data_obj.read_map(1));
        figure
        plot_map(z_input,0,1,data_obj.ROI)
        rectangle = input('Input a fitting rectangle to cut out from the OPM! Format [xUpLeft yUpLeft xDownRight yDownRight]');
        close
        data_info.rectangle=rectangle;
        
        try
            save_info(data_path,data_info)
        catch
            disp(['Saving the rectangle parameter to data_info in ' data_path 'was not possible'])
        end
        
        data_info.RectangleROI = getRectangleROI(data_info.rectangle,data_obj.ROI);
        data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
        
    elseif isfield(data_info, 'rectangle') 
        
        data_info.RectangleROI = getRectangleROI(data_info.rectangle,data_obj.ROI);
        
    else
        data_info.rectangle=[0 0 size(data,1) size(data,2)];
        data_info.RectangleROI = data_obj.ROI;
    end

end

function data = convertKalatskyData(data,data_info)
    data = extractKalatskyDataStruct(data,data_info);
    if length(data_info.stim_order) == 1
        % check if zero
        if data_info.stim_order(1) == 0
            return
        else
            error('Kalatsky data conversion not possible for non-zero single stim!')
        end
    else
        data = extractStimResponse(data,data_info.stim_order);
    end
end

function data_stim = extractStimResponse(data,stim_order)
    data_stim = zeros(size(data,1),size(data,2),length(stim_order),size(data,4));
    for ii = 1:size(data,4)
        z = data(:,:,1,ii);
        for jj = 1:size(stim_order)
            data_stim(:,:,jj,ii) = real(z.*ones(size(z)).*exp(-2i*stim_order(jj)/180*pi));
        end
    end
end

function data_array = extractKalatskyDataStruct(data,data_info,FieldPrefix)
    if nargin < 3
        FieldPrefix = {'cw' 'aw'};
    end

    %size(data_array) = [x,y,1,num struct]
    N_samples = getN_samples(data,FieldPrefix);
    data_array = zeros(data_info.field_size_pix(1),data_info.field_size_pix(2),1,N_samples);

    FieldExists = true;
    ii = 1;
    while FieldExists
        for jj = 1:length(FieldPrefix)
            FieldName = [FieldPrefix{jj} num2str(ii)];
            FieldExists = isfield(data,FieldName);
            if FieldExists
                data_array(:,:,1,ii) = data.(FieldName);
            end
        end
        ii = ii + 1;
    end
end

function BloodVesselImg = getBloodVesselImgFromMean(Data)
    BloodVesselImg = mean(Data,[3,4]);
end

function N_samples = getN_samples(data,FieldPrefix)
    N_samples = 0;
    FieldExists = true;
    ii = 1;
    while FieldExists
        for jj = 1:length(FieldPrefix)
            FieldName = [FieldPrefix{jj} num2str(ii)];
            FieldExists = isfield(data,FieldName);
            if FieldExists
                N_samples = N_samples + 1;
            end
        end
        ii = ii + 1;
    end
end

function ROI = getROI(data_path,ID)
    load([data_path ID '_Mask.mat'],'a')
    ROI = a;
end
