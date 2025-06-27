clear all
close all

%% add path to dir with processing functions
addpath 'OPM_Processing'


%% set data folder
main_folder = 'MinimalExample/';
ID = 'F10-016';
getCI = true; % get confidence intervals
SizeGaussKernelPwDensityCalc = 0.5; % relative size of the kernel for the pinwheel density calculation
Confidence = 0.95; % confidence level for the confidence intervals
scale = 100; % scale for the noise covariance calculation
number_bootstrapsamples = 100; % number of bootstrap samples

%% define result directory
result_dir = [main_folder 'Results/' ID '/'];

%% load data
data_file = [result_dir 'ProcessedData.mat'];
if exist(data_file, 'file')
    load(data_file,'data');
else
    error('Data file not found: %s', data_file);
end
disp('Dimnsions of data:');
disp(size(data));

%% load ROI
data_info_file = [result_dir 'data_info.mat'];
if exist(data_info_file, 'file')
    load(data_info_file,'ROI');
else
    error('ROI file not found: %s', data_info_file);
end
disp('Dimnsions of ROI:');
disp(size(ROI));


%% load data_info
data_info_file = [result_dir 'data_info.mat'];
if exist(data_info_file, 'file')
    load(data_info_file,'data_info');
else
    error('data_info file not found: %s', data_info_file);
end
disp(data_info)

%% exztract parameters from data_info
animal = data_info.animal;
analysis_range_mm = data_info.settings.analysis_range_mm;
lowpass_cutoff_range_mm =data_info.settings.lowpass_cutoffs_mm;
cleaning_method = data_info.settings.data_cleaning_method;
filter_lowpass_cutoff_mm = data_info.settings.lowpass_mm;
filter_highpass_cutoff_mm = data_info.settings.highpass_mm;


%% prepare result directory
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

%% prepare data object
data_obj = data_handle_corrected(data_info,data,ROI);
data_obj.set_data_cleaning_method(cleaning_method);
data_obj.prepare_samples_array(number_bootstrapsamples);
data_obj.set_filter_parameters('lowpass',filter_lowpass_cutoff_mm)
data_obj.set_filter_parameters('highpass',filter_highpass_cutoff_mm)



%% get column spacing filtered
disp('get column spacing filtered')
smallest_w_mm = min(analysis_range_mm);
largest_w_mm = max(analysis_range_mm);
w_step_mm = mean(diff(analysis_range_mm));
[mean_spacing_mm,local_spacing_mm,newROI] = getColumnsSpacing(data_obj, ...
    result_dir,smallest_w_mm,largest_w_mm,w_step_mm,getCI,true);


%% get pinwheel infos
disp('get pinwheel infos')
do_plotting=0;
PwInfo = getPinwheelInfos(data_obj, ...
    local_spacing_mm,result_dir,newROI,getCI,Confidence, ...
    do_plotting,lowpass_cutoff_range_mm,SizeGaussKernelPwDensityCalc);

%% get CI filtered
disp('get CI filtered')
DoFilter = true;
calcCIs(data_obj,Confidence,DoFilter,result_dir);

%% get CI unfiltered
disp('get CI unfiltered')
DoFilter = false;
calcCIs(data_obj,Confidence,DoFilter,result_dir);


%% testModularityOPM
disp('testModularityOPM')
testModularityOPM(data_obj,result_dir,mean_spacing_mm, ...
    analysis_range_mm,number_bootstrapsamples)



%% testPWsOPM
disp('testPWsOPM')
testPWsOPM(data_obj,PwInfo.pinwheel_stats, ...
    number_bootstrapsamples,result_dir)



%% get Noise Covarienaces unfiltered
disp('get Noise Covarienaces unfiltered')
DoFilter = false;
getNoiseCovariances(data_obj,result_dir,'vector',DoFilter,scale);
getNoiseCovariances(data_obj,result_dir,'align',DoFilter,scale);

%% get Noise Covarienaces filtered
disp('get Noise Covarienaces filtered')
DoFilter = true;
getNoiseCovariances(data_obj,result_dir,'vector',DoFilter,scale);
getNoiseCovariances(data_obj,result_dir,'align',DoFilter,scale);



disp(['mean spacing [mm] ' num2str(mean_spacing_mm)])
disp(['mean pw density ' num2str(PwInfo.MeanPwDensity)])
disp(['mean pw number ' num2str(PwInfo.NumberPw)])



disp('Analysis Finished!')

disp('plot results')


%% plot factsheet
FigureFile = [DataFolder 'Factsheet_' animal ' ' data_info.ID];
disp(FigureFile)
if getCI
    PlotFactSheetPage(data_info,data_obj,DataFolder,FigureFile)
else
    PlotResultWoCI(animal,data_info,data_obj,DataFolder,FigureFile)
end