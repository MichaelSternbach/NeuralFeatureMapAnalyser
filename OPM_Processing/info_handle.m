function [this_info, file_path] = info_handle(data_set, set_ID, mount_point_data, data_info_file_list)
    % Default initialisation
    this_info = [];
    file_path = [];

    % Default mount point
    if nargin < 3 || isempty(mount_point_data)
        mount_point_data = '~/CIDBN/';
    end

    % Default CSV file path (same directory as this function)
    if nargin < 4 || isempty(data_info_file_list)
        function_dir = fileparts(mfilename('fullpath')); % Get function directory
        data_info_file_list = fullfile(function_dir, 'experiment_info.csv');
    end

    % Check if the CSV file exists
    if ~isfile(data_info_file_list)
        error('CSV file not found at: %s', data_info_file_list);
    end

    % Read CSV data
    data = readtable(data_info_file_list);

    % Filter data for the selected dataset
    matches = strcmpi(data.Dataset, data_set);
    if ~any(matches)
        error('Dataset "%s" not recognised. Check the CSV file.', data_set);
    end

    % Get dataset information
    selected_data = data(matches, :);
    make_info_rel_path = selected_data.MakeInfo{1};
    destination_folder_rel_path = selected_data.DestinationFolder{1};

    make_info = fullfile(mount_point_data, make_info_rel_path);
    destination_folder = fullfile(mount_point_data, destination_folder_rel_path);

    % Display available experiments if no set_ID is provided
    if nargin < 2 || isempty(set_ID)
        disp('Available experiments:');
        disp(table(selected_data.ExperimentID, selected_data.ExperimentName));
        this_info = table(selected_data.ExperimentID, selected_data.ExperimentName);
        return;
    end

    % Find the specific experiment
    if isnumeric(set_ID)
        experiment_row = selected_data(selected_data.ExperimentID == set_ID, :);
    else
        error('set_ID must be numeric.');
    end

    if isempty(experiment_row)
        error('Experiment ID %d not found for dataset "%s".', set_ID, data_set);
    end

    % Construct file path and load info
    experiment_name = experiment_row.ExperimentName{1};
    file_path = fullfile(destination_folder, experiment_name, '/');
    info_path = fullfile(destination_folder, experiment_name, 'exp_info');

    if isfile([info_path '.mat'])
        disp(['Loading data info from: ', info_path, '.mat']);
        tmp = load([info_path '.mat'], 'info');
        this_info = tmp.info;
    elseif isfile([info_path '.json'])
        disp(['Loading data info from: ', info_path, '.json']);
        this_info = readJSONToStruct([info_path '.json']);
    else
        if isMatlabScript(make_info)
            disp(['Info file not found at: %s', info_path, '.']);
            disp(['run ', make_info]);
            run(make_info);
        elseif isJSON(make_info)
            disp(['Info file not found at: %s', info_path, '.']);
            disp(['Reading JSON file: ', make_info]);
            this_info = readJSONToStruct(make_info);
        else
            error('Unknown info file format: %s', make_info);
        end

    end

    % add info path
    this_info.info_path = info_path;

end

function isJSON = isJSON(file_path)
    % Check if a file is a JSON file
    [~, ~, ext] = fileparts(file_path);
    isJSON = strcmpi(ext, '.json');
end

function is_script = isMatlabScript(file_path)
    % Check if a file is a MATLAB script
    [~, ~, ext] = fileparts(file_path);
    is_script = strcmpi(ext, '.m');
end

function is_csv = isCSV(file_path)
    % Check if a file is a CSV file
    [~, ~, ext] = fileparts(file_path);
    is_csv = strcmpi(ext, '.csv');
end

% function data_info_structs = read_csv_to_struct(csv_file)
%     % Reads a CSV file and outputs the data in the original MATLAB struct format
%     % Input:
%     %   csv_file - the path to the CSV file
%     % Output:
%     %   data_info_structs - array of structs with the data

%     % Read the CSV file
%     opts = detectImportOptions(csv_file);
%     opts = setvaropts(opts, {'stim_order', 'expIds', 'refWin', 'sigWin', 'partId', 'experiment_num'}, 'Type', 'char');
%     data_table = readtable(csv_file, opts);

%     % Initialise an array of structs
%     num_rows = height(data_table);
%     data_info_structs = struct([]);

%     for i = 1:num_rows
%         % Parse field_size_pix
%         field_size_pix = [data_table.field_size_pix_x(i), data_table.field_size_pix_y(i)];
        
%         % Parse stim_order as numeric array
%         stim_order = str2double(split(data_table.stim_order{i}, ';'));
        
%         % Parse expIds as cell array of numeric arrays
%         expIds_str = split(data_table.expIds{i}, ';');
%         expIds = cellfun(@(x) str2num(x), expIds_str, 'UniformOutput', false); %#ok<ST2NM>
        
%         % Parse refWin and sigWin as numeric arrays
%         refWin = eval(data_table.refWin{i}); %#ok<EVLK>
%         sigWin = eval(data_table.sigWin{i}); %#ok<EVLK>
        
%         % Parse partId as char array
%         partId = split(data_table.partId{i}, ';');
        
%         % Parse experiment_num as numeric array
%         experiment_num = eval(data_table.experiment_num{i}); %#ok<EVLK>
        
%         % Create the struct
%         data_info_structs(i).ID = data_table.ID{i};
%         data_info_structs(i).animal = data_table.animal{i};
%         data_info_structs(i).field_size_pix = field_size_pix;
%         data_info_structs(i).pix_per_mm = data_table.pix_per_mm(i);
%         data_info_structs(i).stim_order = stim_order;
%         data_info_structs(i).expIds = expIds;
%         data_info_structs(i).refWin = refWin;
%         data_info_structs(i).sigWin = sigWin;
%         data_info_structs(i).partId = partId;
%         data_info_structs(i).weight_in_grams = data_table.weight_in_grams(i);
%         data_info_structs(i).age_days = data_table.age_days(i);
%         data_info_structs(i).gender = data_table.gender{i};
%         data_info_structs(i).date_recording = data_table.date_recording{i};
%         data_info_structs(i).settings.lowpass_mm = data_table.settings_lowpass_mm(i);
%         data_info_structs(i).settings.highpass_mm = data_table.settings_highpass_mm(i);
%         data_info_structs(i).experiment_num = experiment_num;
%     end
% end


function [this_info,file_path] = info_handle_old(data_set,set_ID,data_dir)
this_info=[];
file_path=[];

% if no arguments, display possible animals
if nargin == 0
    disp('The available data sets are:')
    disp({'ferret';'pairing';'chronic';'cat';'macaque_sam'})
    return
end
if nargin < 3
    data_dir ='~/CIDBN1/';
end
% check which data set is called in animal
switch lower(data_set)
    case {'ferret','ferrets'}
        make_info = strcat(data_dir,'neurodyn/PairingData/Ferret_Whitney/ISI_Data/make_info_ferret.m');
        destination_folder = strcat(data_dir,'neurodyn/PairingData/Ferret_Whitney/Analysis/');
        experiment_IDs = {1,'F09-191';2,'F10-016';3,'F10-017';4,'F10-040';5,'F10-041';...
            6,'F10-042';7,'F10-046';8,'F10-047';9,'F10-048';10,'F10-049';11,'F10-066';...
            12,'F10-076';13,'F10-078';14,'F10-097';15,'F10-117';...
            16,'F10-118';17,'F10-121';18,'F10-122';19,'F10-139';20,'F10-140';21,'F10-157';...
            22,'F10-158';23,'F10-159';24,'F10-160';25,'F10-165';26,'F10-166';27,'F10-190';...
            28,'F10-191';29,'F10-192';30,'F10-193'};
    case {'wallaby'}
        make_info = [data_dir 'WallabyJung/make_info_Wallaby.m'];
        destination_folder = [data_dir 'WallabyJung/'];
        experiment_IDs = {1,'wallabyB';2,'wallabyC';3,'wallabyD';4,'wallabyE';5,'wallabyF';6,'wallabyH'};
        
    case {'dunnart'}
        make_info = [data_dir 'DunnartJung/make_info_Dunnart.m'];
        destination_folder = [data_dir 'DunnartJung/'];
        experiment_IDs = {1,'DunnartAH';2,'DunnartAN';3,'DunnartAN_RightHemisphere';4,'DunnartAM';5,'DunnartAP';6,'dunnartQ';7,'dunnartR';8,'dunnartS';9,'dunnartT';10,'DunnartXX'};
    case {'cat','cats'}
        make_info = strcat(data_dir,'neurodyn/PairingData/Cat_Loewel/Data/make_info_cat.m');
        destination_folder = strcat(data_dir,'neurodyn/PairingData/Cat_Loewel/Analysis/ICMS/');
        experiment_IDs = {1,'K058';2,'K060';3,'K108';4,'K109';5,'K114';6,'K116';...
            7,'K119';8,'K120';9,'K128';10,'K239';11,'K242';12,'K333';13,'K341';14,'K343'};
    case {'cats_all'}
        make_info = '';
        destination_folder = strcat(data_dir,'neurodyn/PairingData/Cat_Loewel/Analysis/Cats/');
        experiment_IDs = {1,'K001';2,'K002';3,'K003';4,'K004';5,'K005';6,'K006';7,'K007';...
            8,'K008';9,'K009';10,'K010';11,'K011';12,'K012';13,'K013';14,'K014';15,'K015';...
            16,'K058';17,'K060';18,'K108';19,'K109';20,'K114';21,'K115';22,'K116';23,'K119';...
            24,'K120';25,'K122';26,'K124';27,'K125';28,'K128';29,'K129';30,'K131';31,'K136';...
            32,'K239';33,'K242';34,'K256';35,'K333';36,'K341';37,'K343';38,'K555'};
    case {'pairing'}
        % info_path = '/pairing/PairingData/Ferret_Whitney/metadata/';
        destination_folder =  '/pairing/PairingData/Ferret/ISI_Data/';
        experiment_IDs = {1,'F10-040';2,'F10-041';3,'F10-042';4,'F10-046';5,'F10-047';...
            6,'F10-048';7,'F10-049';8,'F10-066';9,'F10-076';10,'F10-078';11,'F10-097';...
            12,'F10-118';13,'F10-121';14,'F10-122';15,'F10-139';16,'F10-140';...
            17,'F10-158';18,'F10-159';19,'F10-160';20,'F10-166';21,'F10-190';...
            22,'F10-191';23,'F10-192';24,'F10-193'};
    
    case {'chronic'}
        make_info = '/pairing/PairingData/Ferret_Whitney/ISI_Data/make_info_ferret_chronic.m';
        destination_folder = 'neurodyn/PairingData/Ferret_Whitney/Chronic_analysis/';
        experiment_IDs = {1,'F09-174_LCtx';2,'F10-009_LCtx';3,'F10-163_LCtx';4,'F10-163_RCtx';...
            5,'F10-179_LCtx';6,'F10-179_RCtx';7,'F10-189_LCtx';8,'F10-189_RCtx';...
            9,'F10-198_LCtx';10,'F10-198_RCtx';11,'F10-199_LCtx';12,'F10-199_RCtx'};
    case {'macaque_sam'}
        make_info = strcat(data_dir,'neurodyn/PairingData/Macaque_Angelucci/Data/make_info_macaque.m');
        destination_folder = strcat(data_dir,'neurodyn/PairingData/Macaque_Angelucci/Analysis/');
        experiment_IDs = {1,'MK319LH';2,'MK327LH';3,'MK356RH';4,'MK364LH';5,'MK368RH';6,'MK373LH';7,'MK374RH';8,'MK368LH';9,'MK365LH'};

    case {'macaque_ikezoe'}
         make_info = strcat(data_dir,'neurodyn/PairingData/Macaque_Okamoto/make_info_files.m');
         destination_folder = strcat(data_dir,'neurodyn/PairingData/Macaque_Okamoto/Analysis/');
         experiment_IDs = {1,'monkey1';2,'monkey2';3,'monkey3';4,'monkey4'};

    case {'galago'}
        make_info = strcat(data_dir,'neurodyn/PairingData/Galago_White/make_info_galago.m');
        destination_folder = strcat(data_dir,'neurodyn/PairingData/Galago_White/Analysis/');
        experiment_IDs = {1, 'GC9516';2, 'GC9601';3 ,'GC9602';4 ,'GC9603';...
            5 ,'GC9604';6 ,'GC9605';7 ,'GC9606';8, 'GC9607';9, 'GC9701';...
            10, 'GC9702';11, 'GC9703'};
    
    case {'microcebus','mouse lemur'}
        make_info = strcat(data_dir,'neurodyn/PairingData/Microcebus_Huber/make_info_microcebus.m');
        destination_folder = strcat(data_dir,'neurodyn/PairingData/Microcebus_Huber/Analysis/');
        experiment_IDs = {1,'Chip'; 2,'Dale'; 3,'Argi'; 4,'Burrito'; 5 ,'Hashtag'};
        
    case {'galago_casagrande'}
        make_info = '/home/franzj/num/matlb/Raw_data_analysis/make_info_casagrande.m';
        destination_folder = '/home/franzj/Casagrande/Analysis/';
        experiment_IDs = {1,'PrimateData08.07.01';2,'PrimateData05.28.02';...
            3,'PrimateData07.24.02';4,'PrimateData08.06.02';...
            5,'PrimateData06.10.03';6,'PrimateData06.24.03';7,'PrimateData07.14.03'};
        
    case {'aotus'}
        make_info = '/home/franzj/num/matlb/Raw_data_analysis/make_info_casagrande.m';
        destination_folder = '/home/franzj/Casagrande/Analysis/';
        experiment_IDs = { 1,'PrimateData11.07.01';2,'PrimateData11.01.01';...
            3,'PrimateData08.21.01';4,'PrimateData12.06.01';5,'PrimateData12.12.01';...
            6, 'PrimateData01.15.02';7,'PrimateData02.12.02';8,'PrimateData02.26.02';...
            9,'PrimateData04.16.02';10,'PrimateData04.24.02';11, 'PrimateData06.18.02'};
        
        %%unclassified
        %'PrimateData03.19.02';'PrimateData08.13.03';'PrimateData07.17.01';
        %'PrimateData06.08.01'probably only a testing folder
        
    otherwise
        error('animal entry not recognized')
end

% If only one argument, return a list of availalbe experiments
if nargin==1
    disp('Available experiments:')
    disp(experiment_IDs)
    this_info = experiment_IDs;
    return
end

% return path and info
if isnumeric(set_ID)
    % return path and info of a specific experiment
    %run(make_info)
    file_path = [destination_folder,experiment_IDs{set_ID,2},'/'];
    tmp = load([destination_folder,experiment_IDs{set_ID,2},'/exp_info.mat'],'info');
    %tmp = load([destination_folder,experiment_IDs{set_ID,2},'/info.mat']);
    this_info = tmp.info;
    return   
else
    % make info
    run(make_info)
end

end
