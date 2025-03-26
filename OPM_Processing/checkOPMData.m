function checkOPMData(animals,data_info_file,csv_save_path,AnimalDataFolder)
    
    %% run locally
    % checkOPMData('dunnart:9|wallaby:6','experiment_info.csv','exp_info_fullTest.csv')
    
    %% run on HPC
    % rsync -avvH /home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/OPM_Processing/checkOPMData sternbach1@login-dbn02.hpc.gwdg.de:/scratch/users/sternbach1/OPM_DataPipeline2/
    % sbatch -p cidbn -A cidbn_legacy -c 1 -t 100:00:00 --mem=100G -o /scratch/users/sternbach1/OPM_DataPipeline2/ScriptCheckData_output.txt --export=MODULE="matlab-mcr/R2022b",WORKDIR="/scratch/users/sternbach1/OPM_DataPipeline2/",EXECUTABLE="./checkOPMData",ARGS="dunnart:9|ferret:1 /scratch/users/sternbach1/OPM_DataPipeline2/experiment_info.csv /scratch/users/sternbach1/OPM_DataPipeline2/info_table_Dunnart.csv /home/uni08/cidbn1/" /scratch/users/sternbach1/OPM_DataPipeline2/run_CompiledMatlabCode.sh
    
%     animals.dunnart = 9;
%     data_info_file = '/scratch/users/sternbach1/OPM_DataPipeline2/experiment_info.csv';
%     csv_save_path = '/scratch/users/sternbach1/OPM_DataPipeline2/info_table_Dunnart.csv';
%     AnimalDataFolder = '/home/uni08/cidbn1/';

    if nargin < 4
        AnimalDataFolder = '~/CIDBN/';
    end

    animals = convertChar(animals);

    FullInfoTable = table();
    for animal = fieldnames(animals)
        experiment_num_list = 1:animals.(animal{1});
        for experiment_num = experiment_num_list
            
            %% load experiment data
            [data_info,data_path,data_obj,~,~] = getAnimalData(animal{1},experiment_num,AnimalDataFolder,data_info_file);
            data_info.data_path = data_path;

            try
                [~,~,DirectionData] = getDirectionData(data_info,data_path,0);
                data_info.DirectionData = 1;
            catch
                disp('No Direction Data')
                data_info.DirectionData = 0;
            end
            
            %% get full info table
            FullInfoTable = [FullInfoTable; getFullInfoTable(experiment_num,data_info,data_obj)];

        end
    end

    %% save table as csv
    writetable(FullInfoTable,csv_save_path);
    disp('Data saved!')
    
end

function full_info_table = getFullInfoTable(experiment_num_list,data_info,data_obj)
    %% make table of data that contains
    % Dataset,MakeInfo,DestinationFolder,ExperimentID,ExperimentName
    % field_size_pix, pix_per_mm, number of stimuli, number of trials, weight_in_grams, age_days, gender, date_recording
    % if data does not exist leave empty

    %% collect necessary data variables
    Dataset = string(lower(data_info.animal));
    MakeInfo = string(data_info.info_path);
    DestinationFolder = string(data_info.data_path);
    ExperimentID = experiment_num_list;
    ExperimentName = string(data_info.ID);

    field_size_pix = data_info.field_size_pix;
    pix_per_mm = data_info.pix_per_mm;
    number_of_stimuli = sum(~isnan(data_info.stim_order));
    number_of_trials = size(data_obj.data,4);
    DirectionData = data_info.DirectionData;

   %% collect optional data variables
    weight_in_grams = addOptionalFields(data_info,'weight_in_grams');
    age_days = addOptionalFields(data_info,'age_days');
    gender = addOptionalFields(data_info,'gender',"string");
    date_recording = addOptionalFields(data_info,'date_recording',"string");



    %% create table
    full_info_table = table(Dataset,MakeInfo,DestinationFolder,ExperimentID,ExperimentName,field_size_pix,pix_per_mm,number_of_stimuli,number_of_trials,DirectionData,weight_in_grams,age_days,gender,date_recording);

    
end

function field = addOptionalFields(Structure,field,type_)
    if nargin <3
        type_ = 'none';
    end
    if isfield(Structure,field)
        if type_ == "string"
            field = string(Structure.(field));
        else
            field = Structure.(field);
        end
    else
        field = [];
    end
end