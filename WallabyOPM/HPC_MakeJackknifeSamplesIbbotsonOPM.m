function HPC_MakeJackknifeSamplesIbbotsonOPM(Animal,TrialsToUse1, TrialsToUse2)
    datafolder = '/scratch1/users/sternbach1/MarsupialData/Data/Maps/';
    %datafolder = '~/Cloud/PhD/data/data share/Wallaby data/Maps/WallabyH/';
    
    %WallabyList=['H' 'C'];
    
    %TrialsToUse = TrialsToUse(1):TrialsToUse(2);
    
    disp('parameter')
    disp(Animal)
    disp(num2str(TrialsToUse1))
    disp(num2str(TrialsToUse2))
    disp('--------------------------------------------')
    
    disp('load data')
    %load([datafolder 'Wallaby'  WallabyList(Wallaby) '_Image.mat'],'dimg')
    load([datafolder Animal '_Image.mat'],'dimg')
    dimg_cut = dimg(:,str2num(TrialsToUse1):str2num(TrialsToUse2));
    disp(['size' num2str(size(dimg_cut))])
    disp('--------------------------------------------')

    disp('run main function: makeJackknifeIbbotsonOPM')
    disp('--------------------------------------------')
    JackknifeData = makeJackknifeIbbotsonOPM(dimg_cut);
    disp('--------------------------------------------')

    disp('save data')
    Jackknifefilename = [datafolder 'JackknifeSamples' Animal 'Trials' TrialsToUse1 '-' TrialsToUse2 '.mat'];
    save(Jackknifefilename,'JackknifeData');

    disp('--------------------------------------------')
    disp('Finished!')
end
% 
% module load matlab/R2020b
% cd /scratch1/users/sternbach1/MarsupialData/CompiledCode/HPC_MakeJackknifeSamplesIbbotsonOPM
% sbatch -c 1 -N 1 -A all -C scratch --mem=100G -t 48:00:00 -o "/scratch1/users/sternbach1/MarsupialData/output/output_HPC_MakeJackknifeSamplesWallabyH_Trials10.txt" --wrap='matlab -batch "!./HPC_MakeJackknifeSamplesIbbotsonOPM 'H' 1 10"'
