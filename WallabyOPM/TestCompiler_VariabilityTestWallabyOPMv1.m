% module load matlab/R2020b
% cd /scratch1/users/sternbach1/MarsupialData/marsupial-data/WallabyOPM/TestHPC_VariabilityTestWallabyOPMv1
% sbatch -c 1 -N 1 -A all -C scratch --mem=16G -t 1:00:00 -o "/scratch1/users/sternbach1/MarsupialData/output/output_TestHPC_VariabilityTestWallabyOPM.txt" --wrap='matlab -batch "!./TestHPC_VariabilityTestWallabyOPMv1"'

% module load matlab/R2020b
% cd /scratch1/users/sternbach1/MarsupialData/marsupial-data/WallabyOPM
% sbatch -c 1 -N 1 -A all -C scratch --mem=16G -t 1:00:00 -o "/scratch1/users/sternbach1/MarsupialData/output/output_TestHPC_VariabilityTestWallabyOPM.txt" --wrap='matlab -batch "TestHPC_VariabilityTestWallabyOPMv1"'

%addpath('/scratch1/users/sternbach1/MarsupialData/ScriptJason/')

disp('load data')
load('/home/michael/Cloud/PhD/data/data share/Wallaby data/Maps/WallabyH_Image.mat')
disp('--------------------------------------------')

disp('make small test data set')
dimg_small = dimg(:,1:10);
disp('--------------------------------------------')

disp('run main function: VariabilityTestWallabyOPM')
disp('--------------------------------------------')
[AllMaps,MaXDeltaOPM,ConfidenceIntervalls] = VariabilityTestWallabyOPM(dimg_small,2,1);
disp('--------------------------------------------')

disp('save data')
save('TestSmallWallabyH.mat');
disp('--------------------------------------------')
disp('Finished!')