function plotAllDirectionData(DataFolderHPC,ResultFolderHPC,data_info_file_list)
    
    %% run on HPC
    % rsync -avvH /home/michael/Cloud/Cloud/PhD/MarsupialData/OrientationPrefernceMapProcessing/OPM_Processing/plotAllDirectionData sternbach1@login-dbn02.hpc.gwdg.de:/scratch/users/sternbach1/OPM_DataPipeline2/
    % rsync -avvH /home/michael/Cloud/Cloud/PhD/MarsupialData/OrientationPrefernceMapProcessing/OPM_Processing/experiment_info.csv sternbach1@login-dbn02.hpc.gwdg.de:/scratch/users/sternbach1/OPM_DataPipeline2/
    % sbatch -p cidbn -A cidbn_legacy -c 1 -t 100:00:00 --mem=100G -o /scratch/users/sternbach1/OPM_DataPipeline2/plotAllDirectionData_output.txt --export=MODULE="matlab-mcr/R2022b",WORKDIR="/scratch/users/sternbach1/OPM_DataPipeline2/",EXECUTABLE="./plotAllDirectionData",ARGS="/home/uni08/cidbn1/ /scratch/users/sternbach1/OPM_DataPipeline2/DirectionDataPlots/ /scratch/users/sternbach1/OPM_DataPipeline2/experiment_info.csv" /scratch/users/sternbach1/OPM_DataPipeline2/run_CompiledMatlabCode.sh
    
    
    %% parameter
    GIF_TH=4;

    %% check ResultFolderHPC
    if ~isfolder(ResultFolderHPC)
        mkdir(ResultFolderHPC)
    end

    % %% dunnart
    % disp('plot Dunnart')
    % N_animals = 1:9;
    % Data_Set = 'dunnart';
    % plotCIDirectionPreference(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)
    % plotDirectionPreference(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)


    % %% cat
    % disp('plot Cat')
    % N_animals = 1;
    % Data_Set = 'cat_jung';
    % plotCIDirectionPreference(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)
    % plotDirectionPreference(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)

    % %% wallaby
    % disp('plot Wallaby')
    % N_animals = 1:6;
    % Data_Set = 'wallaby';
    % plotCIDirectionPreference(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)
    % plotDirectionPreference(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)

%     %% list data set with # exp as struct
%     data_set.cat_jung = 2;
%     data_set.dunnart = 9;
%     data_set.wallaby = 6;
% 
%     %% loop over data sets
%     data_set_list = fieldnames(data_set);
%     parfor i = 1:length(data_set_list)
%         Data_Set = data_set_list{i};
%         N_animals = 1:data_set.(Data_Set);
%         disp(['plot ' Data_Set])
%         plotOriDirMaps(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)
%         plotDirectionPreference(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)
%         plotCIDirectionPreference(N_animals,Data_Set,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list)
%     end

    experiments.cat_jung = 1:5;
    experiments.microcebus = 1:4;
    experiments.dunnart = 1:6;
    
    plotDirMaps_DunnartPaper(experiments,DataFolderHPC,ResultFolderHPC,GIF_TH,data_info_file_list,0.01,3.5,0.01,100)
    disp('Finished!')
end