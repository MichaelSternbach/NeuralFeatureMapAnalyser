clear all
close all
tic
%% setup data_info
% required
data_info.ID = 'F10-016';
data_info.animal = 'Ferret';
data_info.pix_per_mm = 79.6875;
data_info.stim_order = [NaN, 0, 45, 90, 135];

% optional
data_info.age_in_days = 35;
data_info.weight_in_grams = 103;
data_info.date = '9Mar2010';
data_info.field_size_pix = [510, 510];
data_info.field_size_mm = [6.4, 6.4];


%% set data folder
main_folder = 'MinimalExample/';
data_folder = [main_folder 'Data/' data_info.ID '/'];
result_dir = [main_folder 'Results/' data_info.ID '/'];

%% load data
data_file = [data_folder 'data.mat'];
if exist(data_file, 'file')
    load(data_file,'data');
else
    error('Data file not found: %s', data_file);
end
disp('Dimnsions of data:');
disp(size(data));

%% load ROI
ROI_file = [data_folder 'ROI.mat'];
if exist(ROI_file, 'file')
    load(ROI_file,'ROI');
else
    error('ROI file not found: %s', ROI_file);
end
disp('Dimnsions of ROI:');
disp(size(ROI));


%% prepare result directory
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

%% reduce the number of kernels for parallel computing to prevent memory issues
delete(gcp('nocreate'))
parpool('local', 2);

%% run processin pipeline
removeNanStimSignal = 0; % 0: if NaN stim signal is laready removed, 1: compares methods to remove NaN stim signal
number_bootstrapsamples = 100; % number of bootstrap samples for the analysis (should be at least 100)
data_info = ProcessDataOPM(data,data_info,ROI,result_dir,removeNanStimSignal,number_bootstrapsamples);
toc