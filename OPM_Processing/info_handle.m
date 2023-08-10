function [this_info,file_path] = info_handle(data_set,set_ID)
this_info=[];
file_path=[];

% if no arguments, display possible animals
if nargin == 0
    disp('The available data sets are:')
    disp({'ferret';'pairing';'chronic';'cat';'macaque_sam'})
    return
end
data_dir ='~/CIDBN1/neurodyn/';
% check which data set is called in animal
switch lower(data_set)
    case {'ferret','ferrets'}
        make_info = strcat(data_dir,'/PairingData/Ferret_Whitney/ISI_Data/make_info_ferret.m');
        destination_folder = strcat(data_dir,'/PairingData/Ferret_Whitney/Analysis/');
        experiment_IDs = {1,'F09-191';2,'F10-016';3,'F10-017';4,'F10-040';5,'F10-041';...
            6,'F10-042';7,'F10-046';8,'F10-047';9,'F10-048';10,'F10-049';11,'F10-066';...
            12,'F10-076';13,'F10-078';14,'F10-097';15,'F10-117';...
            16,'F10-118';17,'F10-121';18,'F10-122';19,'F10-139';20,'F10-140';21,'F10-157';...
            22,'F10-158';23,'F10-159';24,'F10-160';25,'F10-165';26,'F10-166';27,'F10-190';...
            28,'F10-191';29,'F10-192';30,'F10-193'};
    case {'wallaby'}
        make_info = '/home/michael/Cloud/PhD/data/data share/Wallaby data/Maps/make_info_Wallaby.m';
        destination_folder = '/home/michael/Cloud/PhD/data/data share/Wallaby data/Maps/';
        experiment_IDs = {1,'WallabyH';2,'WallabyC'};
        
    case {'dunnart'}
        make_info = '/home/michael/Cloud/PhD/data/data share/Dunnart data/make_info_Dunnart.m';
        destination_folder = '/home/michael/Cloud/PhD/data/data share/Dunnart data/';
        experiment_IDs = {1,'DunnartAH';2,'DunnartAN';3,'DunnartAN_RightHemisphere';4,'DunnartAM';5,'DunnartAP';6,'dunnartQ';7,'dunnartR'};
    case {'cat','cats'}
        make_info = '/pairing/PairingData/Cat_Loewel/Data/make_info_cat.m';
        destination_folder = '/pairing/PairingData/Cat_Loewel/Analysis/ICMS/';
        experiment_IDs = {1,'K058';2,'K060';3,'K108';4,'K109';5,'K114';6,'K116';...
            7,'K119';8,'K120';9,'K128';10,'K239';11,'K242';12,'K333';13,'K341';14,'K343'};
    case {'cats_all'}
        make_info = '';
        destination_folder = '/pairing/PairingData/Cat_Loewel/Analysis/Cats/';
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
        destination_folder = '/pairing/PairingData/Ferret_Whitney/Chronic_analysis/';
        experiment_IDs = {1,'F09-174_LCtx';2,'F10-009_LCtx';3,'F10-163_LCtx';4,'F10-163_RCtx';...
            5,'F10-179_LCtx';6,'F10-179_RCtx';7,'F10-189_LCtx';8,'F10-189_RCtx';...
            9,'F10-198_LCtx';10,'F10-198_RCtx';11,'F10-199_LCtx';12,'F10-199_RCtx'};
    case {'macaque_sam'}
        make_info = strcat(data_dir,'PairingData/Macaque_Angelucci/Data/make_info_macaque.m');
        destination_folder = strcat(data_dir,'PairingData/Macaque_Angelucci/Analysis/');
        experiment_IDs = {1,'MK319LH';2,'MK327LH';3,'MK356RH';4,'MK364LH';5,'MK368RH';6,'MK373LH';7,'MK374RH';8,'MK368LH';9,'MK365LH'};

    case {'macaque_ikezoe'}
         make_info = '/pairing/PairingData/Macaque_Okamoto/make_info_files.m';
         destination_folder = '/pairing/PairingData/Macaque_Okamoto/Analysis/';
         experiment_IDs = {1,'monkey1';2,'monkey2';3,'monkey3';4,'monkey4'};

    case {'galago'}
        make_info = '/pairing/PairingData/Galago_White/make_info_galago.m';
        destination_folder = '/pairing/PairingData/Galago_White/Analysis/';
        experiment_IDs = {1, 'GC9516';2, 'GC9601';3 ,'GC9602';4 ,'GC9603';...
            5 ,'GC9604';6 ,'GC9605';7 ,'GC9606';8, 'GC9607';9, 'GC9701';...
            10, 'GC9702';11, 'GC9703'};
    
    case {'microcebus','mouse lemur'}
        make_info = strcat(data_dir,'PairingData/Microcebus_Huber/make_info_microcebus.m');
        destination_folder = strcat(data_dir,'/PairingData/Microcebus_Huber/Analysis/');
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
