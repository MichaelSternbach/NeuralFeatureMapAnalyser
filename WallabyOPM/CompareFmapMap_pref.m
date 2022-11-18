addpath '/home/michael/Cloud/git/vone/MatlabCode/PatchAnalysis'
addpath '/home/michael/Cloud/PhD/data/data share/Wallaby data/Maps'

load('/home/michael/Cloud/PhD/data/data share/Wallaby data/Maps/fmaps_F.mat');
load('/home/michael/Cloud/PhD/data/data share/Wallaby data/Maps/map_pref_F.mat');

theta = (0:7)/8*pi;

addpath('/home/michael/Cloud/PhD/data/data share/Wallaby data/Script')
[phi, r] = oiCalcORMap(fmaps, theta, 'model');

figure()
plot_map(exp(1j*(phi*2-pi)));

figure()
plot_map(exp(1j*map_pref/180*2*pi));

delta = (map_pref/180*2*pi-phi*2);

disp(mean(delta,'all'))
disp(var(delta,1,'all'))