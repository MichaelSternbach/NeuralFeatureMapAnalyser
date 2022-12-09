datafolder = '/scratch1/users/sternbach1/MarsupialData';

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