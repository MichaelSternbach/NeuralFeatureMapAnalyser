function HPC_MakeBootstrapSamplesWallabyOPM(Wallaby,NBootstrapSamples,TrialsToUse1, TrialsToUse2,seed)
    datafolder = '/scratch1/users/sternbach1/MarsupialData/Data/Maps/';
    %datafolder = '~/Cloud/PhD/data/data share/Wallaby data/Maps/WallabyH/';
    
    %WallabyList=['H' 'C'];
    
    %TrialsToUse = TrialsToUse(1):TrialsToUse(2);
    
    disp('parameter')
    disp(Wallaby)
    disp(num2str(NBootstrapSamples))
    disp(num2str(TrialsToUse1))
    disp(num2str(TrialsToUse2))
    disp(num2str(seed))
    disp('--------------------------------------------')
    
    disp('load data')
    %load([datafolder 'Wallaby'  WallabyList(Wallaby) '_Image.mat'],'dimg')
    load([datafolder 'Wallaby'  Wallaby '_Image.mat'],'dimg')
    dimg_cut = dimg(:,str2num(TrialsToUse1):str2num(TrialsToUse2));
    disp(['size' num2str(size(dimg_cut))])
    disp('--------------------------------------------')

    disp('run main function: VariabilityTestWallabyOPM')
    disp('--------------------------------------------')
    AllMaps = VariabilityTestWallabyOPM(dimg_cut,str2num(NBootstrapSamples),str2num(seed));
    disp('--------------------------------------------')

    disp('save data')
    %save([datafolder 'BootstrapSamples' num2str(NBootstrapSamples) 'Wallaby'  WallabyList(Wallaby) '.mat'],'AllMaps');
    save([datafolder 'BootstrapSamples' num2str(NBootstrapSamples) 'Wallaby'  Wallaby 'Trials' TrialsToUse1 '-' TrialsToUse2 '.mat'],'AllMaps');
    disp('--------------------------------------------')
    
    disp('Check if Jackknife data exist:')
    Jackknifefilename = [datafolder 'JackknifeSamplesWallaby'  Wallaby 'Trials' TrialsToUse1 '-' TrialsToUse2 '.mat'];
    if isfile(Jackknifefilename)
        disp('Jackknife data exists!')
    else
        disp('Jackknife data does not exists!')
        disp('--------------------------------------------')
        disp('Start making Jackknife data')
        JackknifeData = makeJackknifeWallabyOPM(dimg_cut);
        disp('--------------------------------------------')
        disp('save data')
        save(Jackknifefilename,'JackknifeData');
    end
    disp('--------------------------------------------')
    disp('Finished!')
end
% 
% module load matlab/R2020b
% cd /scratch1/users/sternbach1/MarsupialData/CompiledCode/HPC_MakeBootstrapSamplesWallabyOPM
% sbatch -c 1 -N 1 -A all -C scratch --mem=100G -t 48:00:00 -o "/scratch1/users/sternbach1/MarsupialData/output/output_HPC_MakeBootstrapSamplesWallabyTrial10.txt" --wrap='matlab -batch "!./HPC_MakeBootstrapSamplesWallabyOPM 'H' 1000 1 10 5"'
