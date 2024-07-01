addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
addpath '/home/michael/Cloud/PhD/data/data share/Wallaby data/Script/'
addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/AGWolfOPMDataPipeline'
addpath '/home/michael/Cloud//git/vone/MatlabCode/PatchAnalysis'

%load('/home/michael/Cloud/PhD/data/data share/Wallaby data/Maps/WallabyH/WallabyH_Image.mat')
%data = PreProcessDataJason(dimg);

%% Jason
theta = [0:7]*pi/8;
maps = mean(data,4);


[op,r] = oiCalcORMap(real(mean(data,4)), [0:7]*pi/8, 'vector');
map1 = r.*exp(2i*op);
% 
figure
plot_map(map1)

%%% Function

% calculate orientation preference by vector averaging
numConds = size(maps,3);
for idx = 1:numConds,
  x(:,:,idx) = maps(:,:,idx) .* cos(2*theta(idx));
  y(:,:,idx) = maps(:,:,idx) .* sin(2*theta(idx));
end

x = sum(x,3);
y = sum(y,3);

phi = atan2(y,x); % atan2 returns angles between -pi and +pi
phi = phi+(sign(phi) < 0)*(2*pi); % 0 < phi < 2*pi

phi = 0.5*phi;


%     yy=[65 220] ;xx =[135 305];; %Contralateral Right Eye
%       xx=[127 282] ;yy =[190 360];; %Ipsilteral Left Eye


%     phi = phi(xx(1):xx(2),yy(1):yy(2));

r = sqrt(x.^2 + y.^2);
%     r= r(xx(1):xx(2),yy(1):yy(2));


%%Chepe


data_set = 'Wallaby';
experiment_num = 1;


[data_info,data_path] = info_handle(data_set,experiment_num);
ROI = true(size(data,1:2));
save([data_path,'/exp_info.mat'],'ROI','-append')

data_info.pixels_per_mm = data_info.pix_per_mm;
data_info.stim_order = theta/pi*180;(0:7)*pi/8;
%data_info.stim_order = data_info.stim_order(1:8);%*0.5;
data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);

map2=data_obj.read_map;%*exp(-0.7231*pi*1i);

figure
plot_map(map2)

% %%% Function Chepe
% 
% sample_number = 1;
% giveVector = false;
% 
% sampleMat = squeeze(data_obj.samples_array(:,:,sample_number));
% 
%  % average the extracted blocks
% map = zeros(size(data_obj.ROI));
% for condition = 1:data_obj.data_parameters.num_stimuli
%     if isnan(data_obj.data_parameters.stimuli_order(condition))
%         continue
%     end
%     tmp_img = squeeze(mean(data_obj.data(:,:,condition,sampleMat(:,condition)),4));
%     map = map + exp(1i*2*pi/180*real(data_obj.data_parameters.stimuli_order(condition)))*tmp_img;
% end


%% Chepe2

stim_order = data_info.stim_order;
ROI=ones(size(data,1:2));
do_gif = false;

map3 = make_map(maps,stim_order,ROI,do_gif);
map3 = normaliseMap(map3);

figure
plot_map(map3)

%%% Function Chepe2

data_3 = mean(data,4);
[nPixY,nPixX,~,nblocks] = size(data_3);

map = zeros(nPixY,nPixX);
for stim_ii = 1:length(stim_order)
    % skip blanks
    if isnan(real(stim_order(stim_ii)))
        continue
    end
    map = map + data_3(:,:,stim_ii) * exp(1i*2*pi/180*real(stim_order(stim_ii)));
end

map = (map-mean(map(ROI),'all'))/std(map(ROI),1,'all');

% disp(size(mean(map(ROI),'all')))
% disp(size(std(map(ROI),1,'all')))

figure
plot_map(map)

% 
% %% Comparison
% %%% Stimulus
% disp(theta)
% disp(data_info.stim_order)
% %%% mean recordings
% ii_x= 50;
% ii_y = 50;
% condition = 4;
% disp(maps(ii_x,ii_y,condition))
% tmp_img = squeeze(mean(data_obj.data(:,:,condition,sampleMat(:,condition)),4));
% disp(tmp_img(ii_x,ii_y))
% 
% disp(sum(abs(tmp_img-maps(:,:,condition)),'all'))

function map = normaliseMap(map)

    %map = (map-mean(map,'all'))/std(map,1,'all');
    map = (map-abs(mean(map,'all')))/abs(std(map,1,'all'));
end
