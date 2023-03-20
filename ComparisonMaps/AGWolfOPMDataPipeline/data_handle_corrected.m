classdef data_handle_corrected < handle
    %% Class that handles experimental data
    %
    % PROPERTIES:
    % data = experimental data; not changed at any time
    % ROI = region of interest. Initially is the whole frame. A ROI file can be
    %       loaded when the class is constructed and it can be changed with an
    %       input ROI
    % data_parameters = specifications about the data, like number of pixels,
    %                   number of blocks, etc
    % filter_parameters = bandpass parameters and padding parameters
    % samples_array = bootstrap samples shuffling data. The first one is
    %                 always the full data.
    %
    % METHODS:
    % data_handle = constructor. Requires as input the path to the data.
    %               Optionally also the path to the ROI, the filter parameters
    %               or the sample array (in order, leave empty if required)
    % set_ROI = Input a ROI array externally
    % set_filter_parameters = modify filter parameters ('parameterName',parameterValue)
    % prepare_samples_array = create array that contains sample indices (numSamples)
    % order_samples_array = organize array given a new positions eg. ([1 5 2 4 3]).
    %                       The first order MUST start with 1.
    % calculate_filter = create bandpass filter using current parameters
    % read_map = read one sample of the data. Input the desired sample.
    %           Optionally, the output can be given as vector using the ROI.
    % filter_map = apply the class filters on the map. Give map as input.
    % array2vector = convert map input from 2D array to vector
    % vector2array = convert map input from vector to 2D array
    %
    % 16/4/2014 chepe@nld.ds.mpg.de
    %%
    
    properties (SetAccess = private)
        data
        ROI
        info
        data_parameters
        filter_parameters
        
        samples_array
        sample1_mean_samples = false
        exported_use = false
        exported_handle = []
        
        reduced_handle = false
        reduced_map
        reduced_condition
        
        lsm_applied = false
        lsm_V = []
        lsm_C = []
        lsm_B = []
        
        typical_scale = []
        
    end
    
    methods
        
        %% Constructor method
        function obj = data_handle_corrected(info,data_input,ROI_input)
            
            % ===========  READ DATA
            %{
              To read the data there are 3 possibilities:
                1) 4D data array [y,x,stim,block]
                2) a string with the path to a file to lead the data array
                3) a cell with multiple strings to load several files and
                combine them
            %}
                       
            if isnumeric(data_input)
                % if option 1) : check that it is 4D and store
                assert(ndims(data_input)==4, 'Data array whould have 4 dimensions.')
                obj.data = data_input;
                blocks_per_file = size(obj.data,4);
                
            else
                % if option 2) -> convert to option 3)
                if ischar(data_input)
                    data_input = {data_input};
                end
                
                blocks_per_file = zeros(size(data_input));
                % -- extract first
                tmp = load(data_input{1},'data');
                total_data = tmp.data;
                blocks_per_file(1)=size(tmp.data,4);
                % -- extract rest
                for data_num = 2:length(data_input)
                    tmp = load(data_input{data_num},'data');
                    total_data = cat(4,total_data,tmp.data);
                    blocks_per_file(data_num)=size(tmp.data,4);
                end
                % add to class
                obj.data = total_data;
            end
            % specifications about data
            obj.data_parameters.pixels_y = size(obj.data,1);
            obj.data_parameters.pixels_x = size(obj.data,2);
            obj.data_parameters.num_stimuli = size(obj.data,3);
            obj.data_parameters.num_blocks = size(obj.data,4);
            obj.data_parameters.blocks_per_file = blocks_per_file;
            
            % ===========  READ ROI
            %{
               To read ROI there are 2 possibilities:
                1) a 2D array with the same dimensions as data
                2) a string with a path to a file with the ROI
                3) empty: all the frame is defined as inside the ROI
            %}
            
            if nargin>2
                if ~ischar(ROI_input)
                    % option 1) a 2D array:
                    assert(all(size(ROI_input)==[obj.data_parameters.pixels_y,obj.data_parameters.pixels_x]),'ROI dimensions do not agree with data.')
                    obj.ROI = logical(ROI_input);
                else
                    % option 2) a path to a file with the ROI:
                    tmp = load(ROI_input,'ROI');
                    if isfield(tmp,'ROI')                                            
                    assert(all(size(tmp.ROI)==[obj.data_parameters.pixels_y,obj.data_parameters.pixels_x]),'ROI dimensions do not agree with data.')
                    obj.ROI = logical(tmp.ROI);
                    else
                        obj.ROI = true(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x);
                    end
                end
            else
                % option 3) take the whole frame as ROI
                obj.ROI = true(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x);
            end
            
            % ===========  READ INFO
            
            %{
              info is a struct that MUST contain:
                - stim_order
                - either field_size_mm or pixels_per_mm
                
              info CAN contain: (to define filter settings)
                - animal 
                - settings     
       
              if info is a 2x1 cell, the second element is a string that 
              specifies which kind of stimuli list to use (.e.g
              'binocular')
            %}
            
            if iscell(info)
                % select the stim_order to use
                info{1}.stim_order = info{1}.stim_order.(info{2});
                info = info{1};
            else
                if isstruct(info.stim_order)
                    tmp = fieldnames(info.stim_order);
                    info.stim_order = info.stim_order.(tmp{1});
                   disp(['No stimuli type given. Using: ',tmp{1}]) 
                end
            end
            obj.info = info;
            
            % add file path for exporting
            if isnumeric(data_input)
                obj.info.data_path = pwd;
                obj.info.data_name = 'data';
            else                
                [pathstr,name] = fileparts(data_input{1});
                obj.info.data_path = pathstr;
                obj.info.data_name = name;
            end
            
            obj.data_parameters.stimuli_order = info.stim_order;
            if isfield(info,'field_size_mm')
                obj.data_parameters.fieldSize = info.field_size_mm;
                obj.data_parameters.pixels_per_mm = max([obj.data_parameters.pixels_y obj.data_parameters.pixels_x]./info.field_size_mm);
            elseif isfield(info,'pixels_per_mm')
                obj.data_parameters.pixels_per_mm = info.pixels_per_mm;
                obj.data_parameters.fieldSize=[obj.data_parameters.pixels_y obj.data_parameters.pixels_x]./info.pixels_per_mm; % of the largest dimension                
            else
                error('Could not determine the size of the frame in mm.')
            end
            
            % fill filter parameters with predefined values
            if isfield(info,'settings')
                obj.filter_parameters.lowpass = info.settings.lowpass_mm;
                obj.filter_parameters.highpass = info.settings.highpass_mm;
                obj.filter_parameters.rise = 0.05;
            else
                switch lower(info.animal)
                    case {'cat','cats'}
                        obj.filter_parameters.lowpass = 0.3;
                        obj.filter_parameters.highpass = 1.6;
                        obj.filter_parameters.rise = 0.05;
                    case {'ferret','ferrets'}
                        obj.filter_parameters.lowpass = 0.4;
                        obj.filter_parameters.highpass = 1.5;
                        obj.filter_parameters.rise = 0.05;
                    case {'galago','galagos'}
                        obj.filter_parameters.lowpass = 0.2;
                        obj.filter_parameters.highpass = 1.2;
                        obj.filter_parameters.rise = 0.05;
                    case {'shrews','shrew','treeshrews','treeshrew'}
                        obj.filter_parameters.lowpass = 0.2;
                        obj.filter_parameters.highpass = 1.2;
                        obj.filter_parameters.rise = 0.05;
                    case {'macaques','macaque'}
                        obj.filter_parameters.lowpass = 0.4;
                        obj.filter_parameters.highpass = 1.2;
                        obj.filter_parameters.rise = 0.05;
                    case {'mouse lemur'}
                        obj.filter_parameters.lowpass = 0.26;
                        obj.filter_parameters.highpass = 1;
                        obj.filter_parameters.rise = 0.05;
                    otherwise 
                        warning('Unknown species: using ferret filter settings')
                        obj.filter_parameters.lowpass = 0.4;
                        obj.filter_parameters.highpass = 1.5;
                        obj.filter_parameters.rise = 0.05;
                end
            end
            
            % ===========  INITIAL CALCULATIONS
            
            % define padding parameters to have square power of 2 array
            pixels = max(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x);
            
            obj.filter_parameters.padd_size_y = 2^nextpow2(pixels);
            obj.filter_parameters.padd_list_y = ...
                ceil((obj.filter_parameters.padd_size_y-obj.data_parameters.pixels_y)/2)...
                +[1:obj.data_parameters.pixels_y];
            
            obj.filter_parameters.padd_size_x = 2^nextpow2(pixels);
            obj.filter_parameters.padd_list_x = ...
                ceil((obj.filter_parameters.padd_size_x-obj.data_parameters.pixels_x)/2)...
                +[1:obj.data_parameters.pixels_x];
            
            % calculate filter
            obj.calculate_filter;
            
            % get initial sample array ready            
            obj.prepare_samples_array(1);
            
        end
        function set_random_generator(obj,seed)
            % initiate a new stream of random numbers. Usefull in the
            % cluster where rng('shuffle') gives the same result
            if nargin==1
                rng('default')
                disp('Random number generator restarted to default stream.')
            elseif isempty(seed)
                rng('shuffle')
                disp('Random number generator restarted using internal clock.')
            else
                rng(seed)
                disp('Random number generator restarted with given seed.')
            end
            
        end
        
        %% Reduce class by removing the data and having only the mean sample
        function reduce_handle(obj)
            % Save mean map, reduce sample array and delete data. This is
            % done to save memory when only the mean is needed
            % WARNING: Not reversible, can cause errors if the class is
            % changed afterwards
            
            % -- extract mean map and condition
            obj.reduced_map = obj.read_map(1);
            obj.reduced_condition = zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_y,obj.data_parameters.num_stimuli);
            for condition = 1:obj.data_parameters.num_stimuli
                obj.reduced_condition(:,:,condition) = obj.read_condition(1,condition);
            end
            
            % -- reduce samples and delete data
            obj.data = [];
            obj.samples_array = squeeze(obj.samples_array(:,:,1));
            
            % -- label as reduced
            obj.reduced_handle = true;
            disp('Handle has been reduced')
        end
        
        %% Clean data
        % Local Similarity Minimization
        function apply_LSM(obj,do_apply,num_templates,window_size_mm)
            % apply_LSM(obj,do_apply,num_templates,window_size_mm)
            % LSM uses the blocks used in samples_array to calculate templates
            % If no input, change state
             
            if nargin<4
                window_size_mm = 0.5;
            end
            if nargin<3                                
                num_templates = 6;
            end
            if nargin<2
                do_apply = ~obj.lsm_applied;
            end
            
            if do_apply
                %disp('LSM method will be applied to maps.')
                obj.lsm_applied = true;
                % prepare variables
                templates = 1:num_templates;
                filt_w = round(window_size_mm * obj.data_parameters.pixels_per_mm); % < hc size = 0.86
                blanks = squeeze(obj.data(:,:,isnan(obj.data_parameters.stimuli_order),unique(obj.samples_array(:,1,1))));
                % get templates
                [obj.lsm_V,obj.lsm_C,obj.lsm_B] = LSM_get_templates(blanks,filt_w,templates);
            else
                disp('LSM method will NOT be applied to maps.')
                obj.lsm_applied = false;
                obj.lsm_V = [];
                obj.lsm_C = [];
                obj.lsm_B = [];
            end
        end
       
        function remove_blocks(obj,remove_list)
            obj.data(:,:,:,remove_list) = [];
            
            obj.data_parameters.num_blocks = size(obj.data,4);
            obj.prepare_samples_array(size(obj.samples_array,3))
        end
        %% Manipulate ROI and filter parameters externally        
        % make new ROI
        function Im = make_ROI(obj)
            % read a clean map
            z = obj.read_map_gif;
            
            % Average normalized difference maps to get image
            Im = zeros(size(z,1),size(z,2));
           
            % Cardinal orientations
            tmp=real(z);
            tmp=(tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
            Im=Im+tmp;
            
            % Oblique orientations
            tmp=imag(z);
            tmp=(tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
            Im=Im+tmp;
            
            % Set range to 0-1
            Im=(Im-min(Im(:)))/(max(Im(:))-min(Im(:)));
            
            % Add previous ROI to activity image as reference
            p=contourc(double(obj.ROI),[1 1]);
            ind = 1;
            while ind<size(p,2)
                line_pix = p(:,ind+[1:p(2,ind)]);
                ind = ind+p(2,ind)+1;
                
                for pix=1:length(line_pix)
                    Im(line_pix(2,pix),line_pix(1,pix))=1;
                end
            end
            
            % == return image without making the ROI
            if nargout==1
               return 
            end
            
            % == open GUI to select ROI
            
            h = figure;
            clf
            obj.ROI = roipoly(Im);
            close(h)
            
        end
        % set the ROI externally
        function set_ROI(obj,ROI)
            
            % reset ROI
            if isempty(ROI)
                obj.ROI = true(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x);
            % check that the size matches
            elseif size(ROI,1)==obj.data_parameters.pixels_y ...
                    && size(ROI,2)==obj.data_parameters.pixels_x
                obj.ROI = logical(ROI);
            else
                disp('WARNING: Size of input ROI does not match frame size.')
                return
            end
        end
        % set filter parameters externally
        function set_filter_parameters(obj,filter_parameter_name,new_value)
            
            % check that field exists
            if ~isfield(obj.filter_parameters, filter_parameter_name)
                disp('WARNING: Name of parameter to change does not exist.')
                return
            else
                obj.filter_parameters.(filter_parameter_name) = new_value;
                obj.calculate_filter;
            end
        end
        
        %% Create and organize sample array
        % set an external sample array
        function set_samples_array(obj,used_samples)
            obj.samples_array = used_samples;
        end
        % create array that holds indices to build samples
        function prepare_samples_array(obj,number_of_samples,blocks2use)
            
            % start random numbers
            obj.set_random_generator('default') % 'shuffle'
            
            % get available blocks to do sampling
            if nargin==2
                blocks2use = 1:obj.data_parameters.num_blocks;
            end
            
            % initialize empty array
            obj.samples_array=zeros(length(blocks2use),...
                obj.data_parameters.num_stimuli,number_of_samples);
            
            % first combination is using all blocks
            for ii=1:obj.data_parameters.num_stimuli
                obj.samples_array(:,ii,1)=blocks2use;
            end
            
            % the rest is filled by random sampling of blocks
            for ii=2:number_of_samples
                for jj=1:obj.data_parameters.num_stimuli
                    obj.samples_array(:,jj,ii)=blocks2use(ceil(length(blocks2use)*rand(length(blocks2use),1)));
                end
            end
        end
        % create jackknife sampling
        function number_of_samples = prepare_jackknife_samples(obj,blocks2use)
            
            % get available blocks to do sampling
            if nargin==1
                blocks2use = 1:obj.data_parameters.num_blocks;
            end
            number_of_samples = length(blocks2use);
            
            % initialize empty array
            obj.samples_array=zeros(length(blocks2use)-1,...
                obj.data_parameters.num_stimuli,number_of_samples);
            
            % fill by removing one data point from each sample
            for ii=1:number_of_samples
                for jj=1:obj.data_parameters.num_stimuli
                    obj.samples_array(:,jj,ii)=blocks2use(setdiff(1:length(blocks2use),ii));
                end
            end
        end
        % create sub_sampling array, where only a number of blocks are used
        % in each sample
        function prepare_subsampling_array(obj,number_of_samples,blocks2use)
            
            % start random numbers
            obj.set_random_generator('default') % 'shuffle'
            
            % if not given, the sample uses all available blocks
            if nargin==2
                blocks2use = obj.data_parameters.num_blocks;
            end
            
            % initialize empty array
            obj.samples_array=zeros(obj.data_parameters.num_blocks,...
                obj.data_parameters.num_stimuli,number_of_samples);
            
            % first combination is using all blocks
            for ii=1:obj.data_parameters.num_stimuli
                obj.samples_array(:,ii,1)=1:obj.data_parameters.num_blocks;
            end
            
            % the rest is filled by random sampling of blocks
            for ii=2:number_of_samples
                % make new block array with the number of blocks to use
                block_array = ceil(obj.data_parameters.num_blocks*rand(blocks2use,1));
                for jj=1:obj.data_parameters.num_stimuli
                    obj.samples_array(:,jj,ii)=block_array(ceil(length(block_array)*rand(obj.data_parameters.num_blocks,1)));
                end
            end
        end
        % Set: mean of samples 2:end as sample 1
        function do_mean_sample1(obj,val)
            obj.sample1_mean_samples = val;
        end
        % organize indices array
        function order_samples_array(obj,new_order)
            
            % check that the new order has the same number of samples as
            % the original array
            if size(obj.samples_array,3)~=length(new_order)
                disp('WARNING: New samples order should have the same size as the original samples array.')
                return
            else
                % check that in the new_order the first sample is allways
                % kept as the original number 1 (original sample)
                if new_order(1)~=1
                    disp('WARNING: First sample in new samples order should be 1.')
                    return
                else
                    obj.samples_array = obj.samples_array(:,:,new_order);
                    disp('Samples array has been re-ordered.')
                end
            end
        end
        % save bootstrap samples
        function export_samples(obj,file_name)
            % make a file with the bootstrap samples such that they are
            % used without calculating them again. Useful in case of
            % cleaning with LSM. File with be saved with the same name as
            % the input, with a _samples suffix.
            
            if nargin==1
                file_name = [obj.info.data_path,'/',obj.info.data_name,'_samples.mat'];
            end
            
            % stop using exported
            obj.exported_use = false;
            obj.exported_handle = [];
            
            % delete previous file is exists
            if exist(file_name,'file')
                delete(file_name)
            end
            
            % make new file
            data = complex(zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,2),0);
            save(file_name,'data','-v7.3');
            
            % add the new maps
            data_file = matfile(file_name,'Writable',true);
            for ii=1:size(obj.samples_array,3)
                data_file.data(:,:,ii) = obj.read_map(ii);
            end
            
            disp(['Bootstrap samples saved in ',file_name])
        end
        % save jackknife samples
        function export_jackknife(obj,file_name)
            % make a file with the bootstrap samples such that they are
            % used without calculating them again. Useful in case of
            % cleaning with LSM. File with be saved with the same name as
            % the input, with a _samples suffix.
            
            if nargin==1
                file_name = [obj.info.data_path,'/',obj.info.data_name,'_jackknife.mat'];
            end
            
            % stop using exported
            obj.exported_use = false;
            obj.exported_handle = [];
            
            % delete previous file is exists
            if exist(file_name,'file')
                delete(file_name)
            end
            
            % make new file
            data = complex(zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,2),0);
            save(file_name,'data','-v7.3');
            
            % add the new maps
            data_file = matfile(file_name,'Writable',true);
            for ii=1:obj.data_parameters.num_blocks
                data_file.data(:,:,ii) = obj.read_map(-ii);
            end
            
            disp(['Jackknife samples saved in ',file_name])
        end
        % use the exported bootstrap samples
        function use_exported_samples(obj,do_it,file_type)
            % file_type = 'sampes' or 'jackknife'
            if nargin==1
                do_it = true;
            end
            if nargin==2
                file_type = 'samples';                
            end
            file_name = [obj.info.data_path,'/',obj.info.data_name,'_',file_type,'.mat'];
            
            % stop using them
            if ~do_it
                % if data was deleted (deployed version), show error
                if isempty(obj.data)
                    error('Class data was deleted! Can''t return to non-exported samples!')
                end
                obj.exported_use = false;
                obj.exported_handle = [];
            else
                % check that file exists
                if ~exist(file_name,'file')
                    disp('ERROR: no file with exported samples found!')
                    obj.exported_use = false;
                    obj.exported_handle = [];
                    return
                else
                    obj.exported_use = true;
                    % if not deployed, lock to the file
                    % otherwise clear data and load samples
                    if ~isdeployed
                       obj.exported_handle = matfile(file_name,'Writable',false);
                       [~,~,num_samples] = size(obj.exported_handle,'data');
                    else
                       obj.data = []; 
                       obj.exported_handle = load(file_name,'data');
                       [~,~,num_samples] = size(obj.exported_handle.data);
                    end
                    disp(['There are ',num2str(num_samples),' samples in the exported file.'])
                end
            end
        end

        %% READING DATA
        % create map from given sample
        function map = read_map(obj,sample_number,giveVector)
            
            % if no input given, choose sample 1
            if nargin==1
                sample_number = 1;
            end
            
            % if not specified, return 2D and not vector
            if nargin<3
                giveVector = false;
            end
            
            % read from exported or create new sample?
            if obj.exported_use
                map = obj.exported_handle.data(:,:,sample_number);
            else
                
                % if requested sample higher than sample mat, give error
                if sample_number == 0 
                    error('Sample number can not be equal zero.')
                elseif sample_number>size(obj.samples_array,3)
                    disp('WARNING: Requested sample higher than available sample mat! Returning sample 1.')
                    sample_number = 1;
                elseif sample_number < -obj.data_parameters.num_blocks
                    error('Sample to remove for Jackknife should be smaller than the number of samples.')
                end
                
                % check if the handle was reduced
                if obj.reduced_handle
                    map = obj.reduced_map;
                    
                    % if required sample is 1, check if mean of data is used or the
                    % mean of the bootstrap samples
                elseif sample_number == 1 && obj.sample1_mean_samples && size(obj.samples_array,3)>1
                    
                    % make an average of all the samples
                    map = zeros(size(obj.ROI));
                    for sample_count=2:size(obj.samples_array,3)
                        map = map + obj.read_map(sample_count,false)/(size(obj.samples_array,3)-1);
                    end
                else
                    
                    % get samples to use
                    if sample_number>0
                        sampleMat = squeeze(obj.samples_array(:,:,sample_number));
                    else
                        % get jackknife sample
                        sampleMat = squeeze(obj.samples_array(:,:,1));
                        sampleMat(abs(sample_number),:) = [];
                    end
                    
                    % average the extracted blocks
                    map = zeros(size(obj.ROI));
                    for condition = 1:obj.data_parameters.num_stimuli
                        if isnan(obj.data_parameters.stimuli_order(condition))
                            continue
                        end
                        tmp_img = squeeze(mean(obj.data(:,:,condition,sampleMat(:,condition)),4));
                        map = map + exp(1i*2*pi/180*real(obj.data_parameters.stimuli_order(condition)))*tmp_img;
                    end
                    % if data is cleaned with LSM, apply to real and imaginary part
                    if obj.lsm_applied
                        map = LSM_apply(real(map),obj.lsm_V,obj.lsm_C,obj.lsm_B) + 1i*LSM_apply(imag(map),obj.lsm_V,obj.lsm_C,obj.lsm_B);
                    end
                end
            end
            % give vector output if requested
            if giveVector
                map = map(obj.ROI);
            end
            
        end
        % create shuffled map
        function map = read_shuffled_map(obj,sample_number,giveVector)
            
            % if no input given, choose sample 1
            if nargin==1
                sample_number = 1;
            end
            
            % if not specified, return 2D and not vector
            if nargin<3
                giveVector = false;
            end
            
            % read map
            map = obj.read_map(sample_number,false);
            
            % phase shuffle
            map(~obj.ROI) = 0;            
            map = ifft2(abs(fft2(map)).*exp(1i*2*pi*rand(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x)));
            map(~obj.ROI) = 0;
            
            % give vector output if requested
            if giveVector
                map = map(obj.ROI);
            end
            
        end
        % read a clean map calculated with Generalized Indicator Functions
        function map = read_map_gif(obj)
            
            % get blocks to use
            blocks2use = unique(squeeze(obj.samples_array(:,1,1)));
                    
            % compress data if binocular and direction
            stim_unique = unique(mod(real(obj.data_parameters.stimuli_order),180)); 
            stim_unique(isnan(stim_unique)) = [];
            
            repeats = sum(bsxfun(@minus,stim_unique,mod(real(obj.data_parameters.stimuli_order),180)')==0,1);
            if ~all(repeats==repeats(1))
                % if each stimulus is not repeated in each eye
                data_clean = zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,length(stim_unique),length(blocks2use));
                for stim_ii=1:length(stim_unique)
                    data_clean(:,:,stim_ii,:) = mean(obj.data(:,:,real(obj.data_parameters.stimuli_order)==stim_unique(stim_ii),:),3);
                end
            else
                data_clean = zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,length(stim_unique),length(blocks2use)*repeats(1));
                for stim_ii=1:length(stim_unique)
                    data_clean(:,:,stim_ii,:) = reshape(...
                        obj.data(:,:,mod(real(obj.data_parameters.stimuli_order),180)==stim_unique(stim_ii),blocks2use),...
                        [obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,length(blocks2use)*repeats(1)]);
                end
            end
                
            % clean data using GIF for each frame separatedly
            data_clean = GIF(data_clean);
            
            % average blocks
            data_clean = squeeze(mean(data_clean,4));
            
            % make maps for each frame
            map = zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x);
            for stim_ii = 1:length(stim_unique)
                map = map + exp(1i*2*pi/180*real(stim_unique(stim_ii)))*data_clean(:,:,stim_ii);
            end
            
        end
        % combine the filter with and without using ROI
        function map = read_chimera(obj,sample_number)
            
            % if no input given, choose sample 1
            if nargin==1
                sample_number = 1;
            end
            
            % read normal filtered map
            A=obj.filter_map(obj.read_map(sample_number));
            
            % read filtered version without removing ROI
            B=obj.filter_map(obj.read_map(sample_number),false);
            
            % combine maps
            dist = bwdist(~obj.ROI);
            dist(dist>5)=5;
            dist = (dist)/5;
            
            % combine maps
            map = A.*dist + B.*(1-dist);
            %map(~obj.ROI)=B(~obj.ROI);
            
            
        end
        % create activity map from given sample
        function map = read_condition(obj,sample_number,condition,giveVector)
            
            % if no input given, choose sample 1
            if nargin==1
                sample_number = 1;
            end
            
            % if no condition specified, return first
            if nargin<3
                condition = 1;
            end
            
            % if not specified, return 2D and not vector
            if nargin<4
                giveVector = false;
            end
            
            % if requested sample higher than sample mat, give error
            if sample_number == 0
                error('Sample number can not be equal zero.')
            elseif sample_number>size(obj.samples_array,3)
                disp('WARNING: Requested sample higher than available sample mat! Returning sample 1.')
                sample_number = 1;
            elseif sample_number < -obj.data_parameters.num_blocks
                error('Sample to remove for Jackknife should be smaller than the number of samples.')
            end
            
            % check if the handle was reduced
            if obj.reduced_handle
                map = obj.reduced_condition(:,:,condition);
                
            elseif sample_number == 1 && obj.sample1_mean_samples && size(obj.samples_array,3)>1
                
                % make an average of all the samples
                map = zeros(size(obj.ROI));
                for sample_count=2:size(obj.samples_array,3)
                    map = map + obj.read_condition(sample_count,condition,false)/(size(obj.samples_array,3)-1);
                end
            else
                
                % get samples to use
                if sample_number>0
                    sampleMat = squeeze(obj.samples_array(:,condition,sample_number));
                else
                    % get jackknife sample
                    sampleMat = squeeze(obj.samples_array(:,condition,1));
                    sampleMat(abs(sample_number),:) = [];
                end
                
                % average the extracted blocks
                map = mean(squeeze(obj.data(:,:,condition,sampleMat)),3);
            end
            
            % give vector output if requested
            if giveVector
                map = map(obj.ROI);
            end
        end
        % create map using certain blocks in the average
        function map = read_blocks(obj,blocks2use)
            map = zeros(size(obj.ROI));
            for condition = 1:obj.data_parameters.num_stimuli
                if isnan(obj.data_parameters.stimuli_order(condition))
                    continue
                end
                map = map + exp(1i*2*pi/180*obj.data_parameters.stimuli_order(condition))*...
                    squeeze(mean(obj.data(:,:,condition,blocks2use),4));
            end
            map = (map - mean(map(obj.ROI)))/std(map(obj.ROI));
        end
        % read OD map
        function map = read_OD(obj,sample_number,giveVector)
            
            % if no input given, choose sample 1
            if nargin==1
                sample_number = 1;
            end
            % if not specified, return 2D and not vector
            if nargin<3
                giveVector = false;
            end
            
            if ~any(imag(obj.data_parameters.stimuli_order)~=0)
                disp('No monocular stimulation available.')
                map = NaN;
                return
            end
            
            % get current sample
            sampleMat = squeeze(obj.samples_array(:,:,sample_number));
            
            % get conditions for monocular stim
            stim_ipsi = imag(obj.data_parameters.stimuli_order)<0 & ~isnan(obj.data_parameters.stimuli_order);            
            stim_contra = imag(obj.data_parameters.stimuli_order)>0 & ~isnan(obj.data_parameters.stimuli_order);
            
            % get samples and average
            data_OD = zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,2,obj.data_parameters.num_blocks);
            for condition = find(stim_contra)
                data_OD(:,:,1,:) = data_OD(:,:,1,:) + obj.data(:,:,condition,sampleMat(:,condition))/sum(stim_contra);                
            end                        
            for condition = find(stim_ipsi)
                data_OD(:,:,2,:) = data_OD(:,:,2,:) + obj.data(:,:,condition,sampleMat(:,condition))/sum(stim_ipsi);
            end
            data_OD = mean(data_OD,4);            
            
            % if data is cleaned with LSM, apply to real and imaginary part
            if obj.lsm_applied
                data_OD = cat(3,...
                    LSM_apply(data_OD(:,:,1),obj.lsm_V,obj.lsm_C,obj.lsm_B),...
                    LSM_apply(data_OD(:,:,2),obj.lsm_V,obj.lsm_C,obj.lsm_B));
            end
            
            % get ocular dominance map
            map = data_OD(:,:,1) - data_OD(:,:,2);
            
            % give vector output if requested
            if giveVector
                map = map(obj.ROI);
            end
        end
        % read ocular dominance map using GIF
        function map = read_OD_gif(obj,method)
            
            if nargin==1
                method = 'cocktail';
            end
            
            if ~any(imag(obj.data_parameters.stimuli_order)~=0)
                    disp('No monocular stimulation available.')
                    map = NaN;
                    return
            end
            
            % get conditions for monocular stim
            stim_ipsi = imag(obj.data_parameters.stimuli_order)<0 & ~isnan(obj.data_parameters.stimuli_order);
            blank_ipsi = imag(obj.data_parameters.stimuli_order)<0 & isnan(obj.data_parameters.stimuli_order);
            
            stim_contra = imag(obj.data_parameters.stimuli_order)>0 & ~isnan(obj.data_parameters.stimuli_order);
            blank_contra = imag(obj.data_parameters.stimuli_order)>0 & isnan(obj.data_parameters.stimuli_order);
            
            blank_all = imag(obj.data_parameters.stimuli_order)==0 & isnan(obj.data_parameters.stimuli_order);
            
            switch method
                case 'all_frames'                    
                    data_OD = zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,2,obj.data_parameters.num_blocks*length(stim_contra));
                    data_OD(:,:,1,:) = reshape(obj.data(:,:,stim_ipsi,:),obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,obj.data_parameters.num_blocks*sum(stim_ipsi));
                    data_OD(:,:,2,:) = reshape(obj.data(:,:,stim_contra,:),obj.data_parameters.pixels_y,obj.data_parameters.pixels_x,obj.data_parameters.num_blocks*sum(stim_contra));
                    
                case 'cocktail'                    
                    if ~isempty(blank_ipsi)
                        data_OD = cat(3,mean(obj.data(:,:,stim_contra,:),3),mean(obj.data(:,:,stim_ipsi,:),3),mean(obj.data(:,:,blank_contra,:),3),mean(obj.data(:,:,blank_ipsi,:),3));
                    elseif ~isempty(blank_all)
                        data_OD = cat(3,mean(obj.data(:,:,stim_contra,:),3),mean(obj.data(:,:,stim_ipsi,:),3),mean(obj.data(:,:,blank_all,:),3));
                    else
                        data_OD = cat(3,mean(obj.data(:,:,stim_contra,:),3),mean(obj.data(:,:,stim_ipsi,:),3));
                    end
            end
            
            % separate sources
            data_OD = squeeze(mean(GIF(data_OD),4));
            
            % get ocular dominance map
            map = data_OD(:,:,1) - data_OD(:,:,2);
            
        end
        % read direction selectivity
        function map = read_DS(obj)
            
            if ~any(real(obj.data_parameters.stimuli_order)>180)
                disp('No direction stimulation available.')
                map = NaN;
                return
            end
            
            % get conditions for monocular stim
            stim_left = real(obj.data_parameters.stimuli_order)<180 & ~isnan(obj.data_parameters.stimuli_order);
            stim_right = real(obj.data_parameters.stimuli_order)>=180 & ~isnan(obj.data_parameters.stimuli_order);
            
            % Make cocktail blank
            data_clean = cat(3,mean(obj.data(:,:,stim_left,:),3),mean(obj.data(:,:,stim_right,:),3));
            
            % separate sources
            data_clean = squeeze(mean(GIF(data_clean),4));
            
            % get ocular dominance map
            map = data_clean(:,:,1) - data_clean(:,:,2);

        end
        
        %% FILTERING AND NORMALIZING
        % function to calculate the filter using given parameters           
        function calculate_filter(obj)
            
            % calculate distances matrix
            [xx,yy] = meshgrid(-obj.filter_parameters.padd_size_y/2:obj.filter_parameters.padd_size_y/2-1 ,...
                -obj.filter_parameters.padd_size_x/2:obj.filter_parameters.padd_size_x/2-1);
            dist=sqrt(xx.^2+yy.^2);
            
            % size of image in mm including PADD
            padd_in_mm = obj.data_parameters.fieldSize .* [obj.filter_parameters.padd_size_y obj.filter_parameters.padd_size_x] ./ [obj.data_parameters.pixels_y obj.data_parameters.pixels_x];
            
            % HighPass
            if isempty(obj.filter_parameters.highpass)
                obj.filter_parameters.filter_highpass=[];
            else
                hp=max(padd_in_mm)/obj.filter_parameters.highpass;
                rise_hp = obj.filter_parameters.rise*hp;
                obj.filter_parameters.filter_highpass=fftshift( 1./(1+exp((dist-hp)./rise_hp)));
            end
            
            % LowPass
            if isempty(obj.filter_parameters.lowpass)
                obj.filter_parameters.filter_lowpass=[];
            else
                lp=max(padd_in_mm)/obj.filter_parameters.lowpass;
                rise_lp = obj.filter_parameters.rise*lp;
                obj.filter_parameters.filter_lowpass=fftshift(1./(1+exp((dist-lp)./rise_lp)));
            end
            
        end
        
%         function calculate_Gaussianfilter(obj)
%             
%             % calculate distances matrix
%             [xx,yy] = meshgrid(-obj.filter_parameters.padd_size_y/2:obj.filter_parameters.padd_size_y/2-1 ,...
%                 -obj.filter_parameters.padd_size_x/2:obj.filter_parameters.padd_size_x/2-1);
%             dist=sqrt(xx.^2+yy.^2);
%             
%             % size of image in mm including PADD
%             padd_in_mm = obj.data_parameters.fieldSize .* [obj.filter_parameters.padd_size_y obj.filter_parameters.padd_size_x] ./ [obj.data_parameters.pixels_y obj.data_parameters.pixels_x];
%             
%             % HighPass
%             if isempty(obj.filter_parameters.highpass)
%                 obj.filter_parameters.Gaussianfilter_highpass=[];
%             else
%                 hp=max(padd_in_mm)/obj.filter_parameters.highpass;
%                 rise_hp = obj.filter_parameters.rise*hp;
%                 std = obj.filter_parameters.highpass*obj.data_parameters.pixels_per_mm;
%                 obj.filter_parameters.Gaussianfilter_highpass=fftshift(exp(-(dist)^2/(2*std*std)) );
%             end
%             %-(x.*x + y.*y)/(2*std*std)
%             
%             % LowPass
%             if isempty(obj.filter_parameters.lowpass)
%                 obj.filter_parameters.Gaussianfilter_lowpass=[];
%             else
%                 lp=max(padd_in_mm)/obj.filter_parameters.lowpass;
%                 rise_lp = obj.filter_parameters.rise*lp;
%                 std = obj.filter_parameters.lowpass*obj.data_parameters.pixels_per_mm;
%                 obj.filter_parameters.Gaussianfilter_lowpass=fftshift(exp(-(dist)^2/(2*std*std)) );
%             end
%             
%         end  
        
        % filter a given file using the defined parameters
        function map = filter_map(obj,map,cut_ROI)
            
            if nargin<2
               map = obj.read_map(1); 
            end
            
            if nargin<3
                cut_ROI = true;
            end
            
            % If input map is vector, check size and convert to 2D
            vector_flag = false;
            if isvector(map)
                vector_flag = true;
                if length(map)~=sum(obj.ROI(:))
                    disp('WARNING: Map input is a vector and does not match the size of the defined ROI')
                    return
                else
                    tmp = zeros(size(obj.ROI));
                    tmp(obj.ROI) = map;
                    map = tmp;
                end
            end
            
            % DC shift
            %map=(map-mean(map(obj.ROI)));
            
            % padd the map and the ROI
            map_padd=zeros(obj.filter_parameters.padd_size_y,obj.filter_parameters.padd_size_x);
            map_padd(obj.filter_parameters.padd_list_y,obj.filter_parameters.padd_list_x)=map;
            ROI_padd=zeros(obj.filter_parameters.padd_size_y,obj.filter_parameters.padd_size_x);
            ROI_padd(obj.filter_parameters.padd_list_y,obj.filter_parameters.padd_list_x)=double(obj.ROI);
            
            %  Apply Fermi Highpass
            if ~isempty(obj.filter_parameters.highpass)
                if cut_ROI
                    map_padd(ROI_padd==0)=0;
                    map_padd = map_padd -  ...
                        fftshift(ifft2(fft2(fftshift(map_padd)).*obj.filter_parameters.filter_highpass))...
                        ./fftshift(ifft2(fft2(fftshift(ROI_padd)).*obj.filter_parameters.filter_highpass));
                    map_padd(ROI_padd==0)=0;
                else
                    map_padd = map_padd -  ...
                        fftshift(ifft2(fft2(fftshift(map_padd)).*obj.filter_parameters.filter_highpass));
                end
            end
            
            % Apply Fermi Lowpass
            if ~isempty(obj.filter_parameters.lowpass)
                if cut_ROI
                    map_padd(ROI_padd==0)=0;
                    map_padd = fftshift(ifft2(fft2(fftshift(map_padd)).*obj.filter_parameters.filter_lowpass))...
                        ./fftshift(ifft2(fft2(fftshift(ROI_padd)).*obj.filter_parameters.filter_lowpass));
                    map_padd(ROI_padd==0)=0;
                else
                    map_padd = fftshift(ifft2(fft2(fftshift(map_padd)).*obj.filter_parameters.filter_lowpass));
                end
            end
            
            % Unpadd and normalize
            map=map_padd(obj.filter_parameters.padd_list_y,obj.filter_parameters.padd_list_x);
            map=(map-mean(map(obj.ROI)))/std(map(obj.ROI));
            %map=map/std(map(obj.ROI));
            
            % If input was a vector, return a vector
            if vector_flag
                map = map(obj.ROI);
            end
        end
        % center the map and normalize
        function map = normalize_map(obj,map)
            
            if nargin<2
               map = obj.read_map(1); 
            end
            
            % DC shift the map again
            map=(map - mean(map(obj.ROI)))/std(map(obj.ROI));
        end
        
        %% Convert map to vector and back
        % convert map from 2D array to vector
        function map = array2vector(obj,map)
            if isvector(map)
                disp('WARNING: Input map is not a 2D array.')
                return
            else
                if ~isequal(size(map),size(obj.ROI))
                    disp('WARNING: Array does not match the size of the defined ROI')
                    return
                else
                    map = map(obj.ROI);
                end
            end
        end
        % convert map from vector to 2D array
        function map = vector2array(obj,map)
            if ~isvector(map)
                disp('WARNING: Input map is not a vector.')
                return
            else
                if length(map)~=sum(obj.ROI(:))
                    disp('WARNING: Vector does not match the size of the defined ROI')
                    return
                else
                    tmp = zeros(size(obj.ROI));
                    tmp(obj.ROI) = map;
                    map = tmp;
                end
            end
        end
        
        %% Get variance of map        
         function var_map = get_sample_variance(obj,do_unfiltered)
            if nargin==1
                do_unfiltered = false;
            end
            if size(obj.samples_array,3)==1
               disp('There is only mean map in sample array!')
                return 
            end
            
            mean_map = obj.filter_map(obj.read_map(1));
            var_map = zeros(obj.data_parameters.pixels_y,obj.data_parameters.pixels_x);
            
            for ii=2:size(obj.samples_array,3)
                if do_unfiltered
                    c_map = obj.read_map(ii);
                else
                    c_map = obj.filter_map(obj.read_map(ii));
                end
                var_map = var_map + ((c_map-mean_map).*conj(c_map-mean_map))/(size(obj.samples_array,3)-1);
            end
            
         end
        
        %% Wrappers to call functions with class        
        % get the radial profile of the power spectrum to define filters
        function [profile,scale_mm,scale_pix] = power_profile(obj,profile_range,do_plot)
            % wrapper to call power_profile.m 
            
            if nargin==1
                do_plot = true;
                profile_range = [10 50];
            end
            if nargin==2
                do_plot = true;
            end
            
            [profile,scale_mm,scale_pix] = power_profile(obj.read_map_gif,obj.ROI,obj.data_parameters.pixels_per_mm,profile_range,do_plot);
        end
        % get typical scale of the map
        function [average_spacing_mm,local_spacing_mm]  = get_column_spacing(obj,smallest_w,largest_w)
            % wrapper to call get_column_spacing.m from the class
            
            % get parameters
            if nargin==1
                switch lower(obj.info.animal)
                    case {'cat','cats'}
                        smallest_w = 0.4;
                        largest_w = 1.5;
                    case {'ferret','ferrets'}
                        smallest_w = 0.4;
                        largest_w = 1.3;
                    case {'galago','galagos'}
                        smallest_w = 0.3;
                        largest_w = 1.1;
                    case {'shrews','shrew','treeshrews','treeshrew'}
                        smallest_w = 0.3;
                        largest_w = 1.1;
                    case {'macaque','macaques'}
                        smallest_w = 0.3;
                        largest_w = 1.1;
                    case {'mouse lemur'}
                        smallest_w = 0.3;
                        largest_w = 1.1;
                    otherwise
                        warning('Unknown species: Calculating from defined filter settings ')
                        lp = obj.filter_parameters.lowpass;
                        if isempty(lp); lp = 0.3; end
                        smallest_w = ceil(10*(lp+0.01))/10;
                        
                        hp = obj.filter_parameters.highpass;
                        if isempty(hp); hp = 1.5; end
                        largest_w = floor(10*(hp-0.01))/10;
                end
            end
            
            % do the calculation
            w_step = 0.05;
            [average_spacing_mm,local_spacing_mm] = get_column_spacing(obj.filter_map(obj.read_map_gif),obj.ROI,obj.data_parameters.pixels_per_mm,smallest_w,largest_w,w_step);      
            
            % save arguments in object
            obj.typical_scale.average_mm = average_spacing_mm;
            obj.typical_scale.average_pixels = average_spacing_mm * obj.data_parameters.pixels_per_mm;
            
            obj.typical_scale.local_mm = local_spacing_mm;
            obj.typical_scale.local_pixels = local_spacing_mm * obj.data_parameters.pixels_per_mm;
            
        end
        % get layout design statistics
        function design_stats = get_design_stats(obj,lowpass_list_mm)
            
            % if the typical scale was not calculated, do it first
            if isempty(obj.typical_scale)
                disp('Calculating typical scale....')
                [~,local_hc_size_mm] = obj.get_column_spacing();
            else
                local_hc_size_mm = obj.typical_scale.local_mm;
            end
            
            % define the range of filter settings to use            
            if nargin==1
                switch lower(obj.info.animal)
                    case {'cat','cats'}
                        lowpass_list_mm = linspace(0.3,1.19,90);
                    case {'ferret','ferrets'}
                        lowpass_list_mm = linspace(0.1,0.99,90);
                    case {'galago','galagos'}
                        lowpass_list_mm = linspace(0.1, 0.79,70);
                    case {'shrews','shrew','treeshrews','treeshrew'}
                        lowpass_list_mm = linspace(0.1, 0.79,70);
                    case {'macaques','macaque'}
                        lowpass_list_mm = linspace(0.1,0.8,71);
                    case {'mouse lemur'}
                        lowpass_list_mm = linspace(0.05, 0.79,70);
                    otherwise
                        warning('Unknown species: Calculating +/- 0.5mm around tyical scale ')
                        lowpass_range = mean(local_hc_size_mm(obj.ROI))+[-0.5 0.5];
                        lowpass_range(1) = max([lowpass_range(1) 0.05]);
                        lowpass_list_mm = linspace(lowpass_range(1),lowpass_range(2),70);
                end                
            end
            
            % calculate stats
            design_stats = get_design_stats(...
                obj.read_map_gif,...
                obj.ROI,...
                local_hc_size_mm,...
                obj.data_parameters.pixels_per_mm,...
                obj.filter_parameters.highpass,...
                lowpass_list_mm);
        end

        %% Plot map
        function [curh,im2plot] = plot_map(obj,map)
            
            %make 2d image
            if size(map,2)==1 && numel(map)>1
                im2d=zeros(size(obj.ROI,1),size(obj.ROI,2));
                if size(map,1)==size(obj.ROI,1)*size(obj.ROI,2)
                    im2d(:)=map;
                else
                    im2d(obj.ROI)=map;
                end
                map=im2d;
            end
            % standard colormaps
            cm = [ 0.2614         0    1.0000;...
                0.4967         0    1.0000;...
                0.7059         0    1.0000;...
                0.9020         0    0.6523;...
                1.0000         0    0.3190;...
                1.0000         0         0;...
                1.0000    0.2353         0;...
                1.0000    0.4314         0;...
                1.0000    0.5882         0;...
                1.0000    0.7451         0;...
                1.0000    0.8824         0;...
                1.0000    1.0000         0;...
                0.8824    1.0000         0;...
                0.7320    1.0000         0;...
                0.5490    1.0000         0;...
                0.2876    1.0000         0;...
                0.1046    1.0000         0;...
                0    1.0000         0;...
                0    0.8562    0.1307;...
                0    0.7059    0.3268;...
                0    0.5490    0.5882;...
                0    0.3922    0.8627;...
                0    0.2092    1.0000;...
                0         0    1.0000];
            mm = [0.3582    0.3582    0.3582;...
                0.3692    0.3692    0.3692;...
                0.3821    0.3821    0.3821;...
                0.3969    0.3969    0.3969;...
                0.4139    0.4139    0.4139;...
                0.4333    0.4333    0.4333;...
                0.4550    0.4550    0.4550;...
                0.4792    0.4792    0.4792;...
                0.5058    0.5058    0.5058;...
                0.5345    0.5345    0.5345;...
                0.5651    0.5651    0.5651;...
                0.5971    0.5971    0.5971;...
                0.6300    0.6300    0.6300;...
                0.6633    0.6633    0.6633;...
                0.6964    0.6964    0.6964;...
                0.7286    0.7286    0.7286;...
                0.7595    0.7595    0.7595;...
                0.7886    0.7886    0.7886;...
                0.8156    0.8156    0.8156;...
                0.8403    0.8403    0.8403;...
                0.8626    0.8626    0.8626;...
                0.8824    0.8824    0.8824;...
                0.8999    0.8999    0.8999;...
                0.9151    0.9151    0.9151;...
                0.9284    0.9284    0.9284;...
                0.9397    0.9397    0.9397;...
                0.9494    0.9494    0.9494;...
                0.9577    0.9577    0.9577;...
                0.9646    0.9646    0.9646;...
                0.9705    0.9705    0.9705;...
                0.9754    0.9754    0.9754;...
                0.9796    0.9796    0.9796;...
                0.9830    0.9830    0.9830;...
                0.9859    0.9859    0.9859;...
                0.9883    0.9883    0.9883;...
                0.9903    0.9903    0.9903;...
                0.9920    0.9920    0.9920;...
                0.9934    0.9934    0.9934;...
                0.9945    0.9945    0.9945;...
                0.9954    0.9954    0.9954;...
                0.9962    0.9962    0.9962;...
                0.9969    0.9969    0.9969;...
                0.9974    0.9974    0.9974;...
                0.9979    0.9979    0.9979;...
                0.9982    0.9982    0.9982;...
                0.9985    0.9985    0.9985;...
                0.9988    0.9988    0.9988;...
                0.9990    0.9990    0.9990;...
                0.9992    0.9992    0.9992;...
                0.9993    0.9993    0.9993;...
                0.9994    0.9994    0.9994;...
                0.9995    0.9995    0.9995;...
                0.9996    0.9996    0.9996;...
                0.9997    0.9997    0.9997;...
                0.9997    0.9997    0.9997;...
                0.9998    0.9998    0.9998;...
                0.9998    0.9998    0.9998;...
                0.9999    0.9999    0.9999;...
                0.9999    0.9999    0.9999;...
                0.9999    0.9999    0.9999;...
                0.9999    0.9999    0.9999;...
                0.9999    0.9999    0.9999;...
                0.9999    0.9999    0.9999;...
                1.0000    1.0000    1.0000];
            
            
            %% angle color
            ori=(angle(map(:)))/(2*pi);ori(ori<0)=1+ori(ori<0);
            
            % find to which color interval it belongs
            intervals=linspace(0,1,2*size(cm,1)+1);
            [~,ori] = histc(ori,intervals(2:2:end));
            ori=ori+1;
            
            % match colors
            anglePlot = (cm(ori(:),:));
            
            %% magnitude color
            
            % scale the magnitude
            magnitude=abs(map);
            max_val = mean(magnitude(obj.ROI))+1*std(magnitude(obj.ROI));
            magnitude=magnitude/max_val; %set to be maximum 1
            magnitude(magnitude>1)=1;
            
            sel=ceil(size(mm,1)*magnitude(:));
            sel(sel==0)=1;
            
            % match color
            magnitudePlot = (mm(sel(:),:));
            
            %% plot
            im2plot = 0.9.*reshape(anglePlot.*magnitudePlot,[size(obj.ROI,1) size(obj.ROI,2) 3]); %repmat(obj.ROI,[1 1 3]).* .*magnitudePlot
            if nargout == 2
                curh = [];
                return
            end
            curh = image(im2plot);
            set(gca,'ydir','reverse')
            
        end
    end
    
end