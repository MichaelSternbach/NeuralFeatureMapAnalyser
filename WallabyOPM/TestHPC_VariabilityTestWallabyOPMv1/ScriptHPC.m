% module load matlab/R2020b
% cd /scratch1/users/sternbach1/MarsupialData/marsupial-data/WallabyOPM/TestHPC_VariabilityTestWallabyOPMv1
% sbatch -c 1 -N 1 -A all -C scratch --mem=16G -t 1:00:00 -o "/scratch1/users/sternbach1/MarsupialData/output/output_TestHPC_VariabilityTestWallabyOPM.txt" --wrap='matlab -batch "ScriptHPC"'


addpath('/scratch1/users/sternbach1/MarsupialData/ScriptJason/')
!./TestHPC_VariabilityTestWallabyOPMv1
