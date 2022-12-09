% module load matlab/R2020b
% cd /scratch1/users/sternbach1/MarsupialData/CompiledCode/v1TestHPC_VariabilityTestWallabyOPM
% sbatch -c 1 -N 1 -A all -C scratch --mem=16G -t 1:00:00 -o "/scratch1/users/sternbach1/MarsupialData/output/output_TestHPC_VariabilityTestWallabyOPM.txt" --wrap='matlab -batch "!./v1TestHPC_VariabilityTestWallabyOPM"'

datafolder = '/scratch1/users/sternbach1/MarsupialData/Data';

disp('load data')
load([datafolder '/Maps/WallabyH_Image.mat'])
disp('--------------------------------------------')

disp('make small test data set')
dimg_small = dimg(:,1:10);
disp('--------------------------------------------')

disp('run main function: VariabilityTestWallabyOPM')
disp('--------------------------------------------')
[AllMaps,MaXDeltaOPM,ConfidenceIntervalls] = VariabilityTestWallabyOPM(dimg_small,2,1);
disp('--------------------------------------------')

disp('save data')
save([datafolder '/Evaluation/VariabilityTest/TestSmallWallabyH.mat']');
disp('--------------------------------------------')
disp('Finished!')