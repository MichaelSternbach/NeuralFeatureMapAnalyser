addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/WallabyOPM'
addpath '/home/michael/Cloud/PhD/data/data share/Wallaby data/Script/'
addpath '/home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/AGWolfOPMDataPipeline'
addpath '/home/michael/Cloud//git/vone/MatlabCode/PatchAnalysis'

load('/home/michael/Cloud/PhD/data/data share/Cat Data/CatBC_Image.mat')
%% pre-process data

refWin = [1:10]; % average over pre-stimulus images (stimulus onset at 2s)

% signal window
sigWin = [31:35]; % FIXME: average around peak of mean intensity plot?
sigWin_ = sigWin;
for i = 1:2

    
%     Combine opposite directions - now looking at 'orientation'
    [aimg(i:2:8,1)] = oiAve([dimg(1:4,:), dimg(5:8,:)]);
     %aimg_(i,1) = oiAve(dimg(9,:));
     
%%     Looking at direction
%     [aimg(i:2:16,1)] = oiAve(dimg(1:8,:));
%     aimg_(i,1) = oiAve(dimg(9,:));
%     
end

%% Cocktail party applied to aimg
aimg_sum =(aimg{1,1}+aimg{2,1}+aimg{3,1}+aimg{4,1}+aimg{5,1}+aimg{6,1}+aimg{7,1}+aimg{8,1})./8;
data = zeros([size(aimg{1,1},1:2) size(aimg,1) size(aimg{1,1},3)]);
    for j = 1:size(aimg,1)
        caimg{j,1}= aimg{j,1}-aimg_sum;
          %%Sign?
        data(:,:,j,:) = -reshape(caimg{j,1},[size(aimg{1,1},1:2) 1 size(aimg{1,1},3)]);
        
    end

data = data(:,:,:,sigWin);    

stim_condition = [0:7]*pi/8;

[op,r] = oiCalcORMap(real(mean(data,4)), stim_condition, 'vector');
map1 = r.*exp(2i*op);

figure
plot_map(map1)

%% Test compatibility with Chepes data class

data_set = 'Wallaby';
experiment_num = 1;


[data_info,data_path] = info_handle(data_set,experiment_num);
ROI = true(size(data,1:2));
save([data_path,'/exp_info.mat'],'ROI','-append')
data_info.pixels_per_mm = data_info.pix_per_mm;
data_info.stim_order = stim_condition/pi*180;
data_info.stim_order = data_info.stim_order(1:8);%*0.5;

data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
map2=data_obj.read_map;%*exp(-0.7231*pi*1i);

figure
plot_map(map2)


%% determine filter settings 1st: manual
power_profile = define_filter_settings(data_info,data_path,data,data_info.stim_order);
% load([data_path,'/exp_info.mat'],'ROI')
% %map = make_map(data,data_info.stim_order,ROI,true);
% %map = make_map(data,data_info.stim_order,ones(size(data,1:2)),false);
% map =data_obj.read_map;

% % figure
% % plot_map(map)
% 
% data_obj = data_handle_corrected(data_info,data,[data_path,'exp_info.mat']);
% disp(data_obj.filter_parameters)
% 
%%% Parameter Chepe
highpass_mm = 3;
lowpass_mm = 1;

% %%% Parameter Jason
% highpass_mm = 0.5;
% lowpass_mm = 0.05;

data_obj.set_filter_parameters('lowpass',lowpass_mm)
data_obj.set_filter_parameters('highpass',highpass_mm)
% %%data_obj.apply_LSM
% 
disp(data_obj.filter_parameters.highpass*data_info.pixels_per_mm)
disp(data_obj.filter_parameters.lowpass*data_info.pixels_per_mm)

z_filtered =data_obj.filter_map(data_obj.read_map);

figure
plot_map(z_filtered)
title('Fermi high- and lowpass Filtered Map');
% 

% %%% recheck power profile
% profile_scale_pixels = 2:2:100;
% profile = calculate_power_profile(map,profile_scale_pixels);
% figure
% plot(profile_scale_pixels,profile)
% figure 
% plot(profile_scale_pixels(1:20),profile(1:20))


% %% determine hyper-column size
% 
% smallest_w = 0.2;
% largest_w = 2; 
% 
% map = data_obj.read_map;
% map_filtered = data_obj.filter_map(map);
% 
% [average_spacing_mm,local_spacing_mm]  = data_obj.get_column_spacing(smallest_w,largest_w);
% disp(average_spacing_mm )

%%% check filter and filtering

% padd the map and the ROI
% map_padd=zeros(data_obj.filter_parameters.padd_size_y,data_obj.filter_parameters.padd_size_x);
% map_padd(data_obj.filter_parameters.padd_list_y,data_obj.filter_parameters.padd_list_x)=map;
% ROI_padd=zeros(data_obj.filter_parameters.padd_size_y,data_obj.filter_parameters.padd_size_x);
% ROI_padd(data_obj.filter_parameters.padd_list_y,data_obj.filter_parameters.padd_list_x)=double(data_obj.ROI);
% 
% %  Apply Fermi Highpass
% if ~isempty(data_obj.filter_parameters.highpass)
%     if cut_ROI
%         map_padd(ROI_padd==0)=0;
%         map_padd = map_padd -  ...
%             fftshift(ifft2(fft2(fftshift(map_padd)).*data_obj.filter_parameters.filter_highpass))...
%             ./fftshift(ifft2(fft2(fftshift(ROI_padd)).*data_obj.filter_parameters.filter_highpass));
%         map_padd(ROI_padd==0)=0;
%     else
%         map_padd = map_padd -  ...
%             fftshift(ifft2(fft2(fftshift(map_padd)).*data_obj.filter_parameters.filter_highpass));
%     end
% end
% 
% % Apply Fermi Lowpass
% if ~isempty(data_obj.filter_parameters.lowpass)
%     if cut_ROI
%         map_padd(ROI_padd==0)=0;
%         map_padd = fftshift(ifft2(fft2(fftshift(map_padd)).*data_obj.filter_parameters.filter_lowpass))...
%             ./fftshift(ifft2(fft2(fftshift(ROI_padd)).*data_obj.filter_parameters.filter_lowpass));
%         map_padd(ROI_padd==0)=0;
%     else
%         map_padd = fftshift(ifft2(fft2(fftshift(map_padd)).*data_obj.filter_parameters.filter_lowpass));
%     end
% end



%%%% Chepe

% figure
% plot_map(map_padd)
% 
% figure
% plot_map(fft2(fftshift(map_padd)))
% 
% figure
% imagesc(fftshift(data_obj.filter_parameters.filter_lowpass)/max(data_obj.filter_parameters.filter_lowpass,[],'all'))
% %%colormap gray;



% %%%% Jason
sigma = 2;
h = fspecial('gaussian', 6*sigma, sigma);
% 
% figure
% imagesc(h/max(h,[],'all'))
% 
% sigma = 20;
% h = fspecial('gaussian', 6*sigma, sigma);
% 
% figure
% imagesc(h/max(h,[],'all'))

%% Comapere Jason and Chepe

%data_obj.calculate_Gaussianfilter()

map = data_obj.read_map;
map_padd=zeros(data_obj.filter_parameters.padd_size_y,data_obj.filter_parameters.padd_size_x);
map_padd(data_obj.filter_parameters.padd_list_y,data_obj.filter_parameters.padd_list_x)=map;

PowerSpectrumMap = abs(fft2(map_padd));
% 
% figure
% imagesc(fftshift(PowerSpectrumMap/max(PowerSpectrumMap,[],'all')))

%disp(sum(data_obj.filter_parameters.filter_lowpass,'all'))

siz = size(map_padd,1);
scale = siz./(1:siz)./data_obj.data_parameters.pixels_per_mm;

% figure
% 
% plot(scale,PowerSpectrumMap(1,:)/max(PowerSpectrumMap,[],'all'),':','DisplayName','power spectrum map')
% xlim([0 max(scale)])
% xlabel('Scale in mm')
% ylabel('Power')

figure

%plot([average_spacing_mm average_spacing_mm],[0 1],'DisplayName','average hypercolumn size')

hold on
plot(power_profile.scale_mm,power_profile.values/max(power_profile.values),':','DisplayName','Power spectrum','LineWidth',2)

hold on

plot(scale,data_obj.filter_parameters.filter_lowpass(1,:)/max(data_obj.filter_parameters.filter_lowpass,[],'all'),'--','DisplayName',['Fermi lowpass ' num2str(data_obj.filter_parameters.lowpass) 'mm'],'LineWidth',1.5)


hold on

plot(scale,data_obj.filter_parameters.filter_highpass(1,:)/max(data_obj.filter_parameters.filter_highpass,[],'all'),'--','DisplayName',['Fermi highpass ' num2str(data_obj.filter_parameters.highpass) 'mm'],'LineWidth',1.5)




%figure
plotFilterJason(2,siz,data_obj.data_parameters.pixels_per_mm,'1st lowpass filter')
hold on
plotFilterJason(20,siz,data_obj.data_parameters.pixels_per_mm,'highpass filter')
hold on
plotFilterJason(8,siz,data_obj.data_parameters.pixels_per_mm,'2nd lowpass filter')

xlabel('Scale in mm')
ylabel('Power')

xlim([0 4])
legend('location','southeast')




% %% determine hyper-column size
% 
% smallest_w = 0.2;
% largest_w = 2;
% 
% 
% 
% 
% 
% 
% 
% map = data_obj.read_map;
% map_filtered = data_obj.filter_map(map);
% 
% [average_spacing_mm,local_spacing_mm]  = data_obj.get_column_spacing(smallest_w,largest_w);
% disp(average_spacing_mm )
% 
% 
% %% determine lowpass filter settings 2nd: pinwheel plateau
% 
% local_w = data_info.pix_per_mm.*local_spacing_mm;
% average_w = data_info.pix_per_mm.*average_spacing_mm;
% 
% z_filtered =data_obj.filter_map(data_obj.read_map);
% 
% figure
% plot_map(z_filtered)
% hold
% border = 10;
% plot([border,border+average_w],[border,border],'color','white')
% 
% pw_density_global = get_pinwheel_density(z_filtered,local_w.*data_obj.ROI);
% 
% lowpass_cutoffs = 0.15:0.01:1.;
% 
% filters = find_lowpass(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs,average_w,local_w);
% 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% High pass filter average image sequence
for i = 1:8
    for k =1:50
        sigma = 20;
        h = fspecial('gaussian', 6*sigma, sigma);
        temp1 = imfilter(caimg{i,1}(:,:,k),h,'symmetric','conv');  
        faimg{i,1}(:,:,k) =caimg{i,1}(:,:,k)-temp1; 
    end
end

%% Low pass filter average image sequence
DataFilteredJason = zeros([size(aimg{1,1},1:2) size(aimg,1) size(sigWin,2)]);
for i = 1:8
    for k =1:50
        sigma = 2; 
        h = fspecial('gaussian', 6*sigma, sigma);
       faimg{i,1}(:,:,k) = imfilter(faimg{i,1}(:,:,k),h,'symmetric','conv');  
    end
    DataFilteredJason(:,:,i,:) = reshape(faimg{i,1}(:,:,sigWin),[size(aimg{1,1},1:2) 1 size(sigWin,2)]);
end

[opFilter,rFilter] = oiCalcORMap(real(mean(DataFilteredJason,4)), stim_condition, 'vector');
mapFilterJason = rFilter.*exp(2i*opFilter);

figure
plot_map(mapFilterJason)
title('Gaussian high- and 1st lowpass filtered map');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% All Gaussian filters


%% High pass filter average image sequence
for i = 1:8
    for k =1:50
        sigma = 20;
        h = fspecial('gaussian', 6*sigma, sigma);
        temp1 = imfilter(caimg{i,1}(:,:,k),h,'symmetric','conv');  
        faimg{i,1}(:,:,k) =caimg{i,1}(:,:,k)-temp1; 
    end
end


%% 1st Low pass filter average image sequence
DataFilteredJason = zeros([size(aimg{1,1},1:2) size(aimg,1) size(sigWin,2)]);
for i = 1:8
    for k =1:50
        sigma = 2; 
        h = fspecial('gaussian', 6*sigma, sigma);
       faimg{i,1}(:,:,k) = imfilter(faimg{i,1}(:,:,k),h,'symmetric','conv');  
    end
    DataFilteredJason(:,:,i,:) = reshape(faimg{i,1}(:,:,sigWin),[size(aimg{1,1},1:2) 1 size(sigWin,2)]);
end


%% 2nd Low pass filter average image sequence
DataFilteredJason = zeros([size(aimg{1,1},1:2) size(aimg,1) size(sigWin,2)]);
for i = 1:8
    for k =1:50
        sigma = 2; 
        h = fspecial('gaussian', 6*sigma, sigma);
       faimg{i,1}(:,:,k) = imfilter(faimg{i,1}(:,:,k),h,'symmetric','conv');  
    end
    DataFilteredJason(:,:,i,:) = reshape(faimg{i,1}(:,:,sigWin),[size(aimg{1,1},1:2) 1 size(sigWin,2)]);
end

[opFilter,rFilter] = oiCalcORMap(real(mean(DataFilteredJason,4)), stim_condition, 'vector');
mapFilterJason = rFilter.*exp(2i*opFilter);


%% Low pass filter average image sequence
DataFilteredJason = zeros([size(aimg{1,1},1:2) size(aimg,1) size(sigWin,2)]);
for i = 1:8
    for k =1:50
        sigma = 8; 
        h = fspecial('gaussian', 6*sigma, sigma);
       faimg{i,1}(:,:,k) = imfilter(faimg{i,1}(:,:,k),h,'symmetric','conv');  
    end
    DataFilteredJason(:,:,i,:) = reshape(faimg{i,1}(:,:,sigWin),[size(aimg{1,1},1:2) 1 size(sigWin,2)]);
end

[opFilter,rFilter] = oiCalcORMap(real(mean(DataFilteredJason,4)), stim_condition, 'vector');
mapFilterJason = rFilter.*exp(2i*opFilter);

figure
plot_map(mapFilterJason)
title('Gaussian high- and 1st and 2nd lowpass filtered map');


% %% Change filter size
% %% Low pass filter average image sequence
% siz=200;
% for i = 1:8
%     for k =1:50
%         sigma = 20;
%         h = fspecial('gaussian', 6*siz, sigma);
%         temp1 = imfilter(caimg{i,1}(:,:,k),h,'symmetric','conv');  
%         faimg{i,1}(:,:,k) =caimg{i,1}(:,:,k)-temp1; 
%     end
% end
% 
% %% Low pass filter average image sequence
% DataFilteredJason = zeros([size(aimg{1,1},1:2) size(aimg,1) size(sigWin,2)]);
% for i = 1:8
%     for k =1:50
%         sigma = 2; 
%         h = fspecial('gaussian', siz, sigma);
%        faimg{i,1}(:,:,k) = imfilter(faimg{i,1}(:,:,k),h,'symmetric','conv');  
%     end
%     DataFilteredJason(:,:,i,:) = reshape(faimg{i,1}(:,:,sigWin),[size(aimg{1,1},1:2) 1 size(sigWin,2)]);
% end
% 
% [opFilter,rFilter] = oiCalcORMap(real(mean(DataFilteredJason,4)), stim_condition, 'vector');
% mapFilterJason = rFilter.*exp(2i*opFilter);
% 
% figure
% plot_map(mapFilterJason)
% title('Filtered Map Jason bigger');





% 
% %% Use chepes filter on Matlab Filter algorithm
% 
% %%% High pass filter average image sequence
% for i = 1:8
%     for k =1:50
%         sigma = 20;
%         hChepe = fftshift(data_obj.filter_parameters.filter_highpass/sum(data_obj.filter_parameters.filter_highpass));%fspecial('gaussian', 6*sigma, sigma);
%         temp1 = imfilter(caimg{i,1}(:,:,k),hChepe,'symmetric','conv');  
%         faimg{i,1}(:,:,k) =caimg{i,1}(:,:,k)-temp1; 
%     end
% end
% 
% %%% Low pass filter average image sequence
% DataFilteredJason = zeros([size(aimg{1,1},1:2) size(aimg,1) size(sigWin,2)]);
% for i = 1:8
%     for k =1:50
%         sigma = 2; 
%         hChepe = fftshift(data_obj.filter_parameters.filter_lowpass/sum(data_obj.filter_parameters.filter_lowpass));%h = fspecial('gaussian', 6*sigma, sigma);
%        faimg{i,1}(:,:,k) = imfilter(faimg{i,1}(:,:,k),hChepe,'symmetric','conv');  
%     end
%     DataFilteredJason(:,:,i,:) = reshape(faimg{i,1}(:,:,sigWin),[size(aimg{1,1},1:2) 1 size(sigWin,2)]);
% end
% 
% [opFilter,rFilter] = oiCalcORMap(real(mean(DataFilteredJason,4)), stim_condition, 'vector');
% mapFilterJason = rFilter.*exp(2i*opFilter);
% 
% figure
% plot_map(mapFilterJason)
% title('Filter Chepe Function Matlab');


function plotFilterJason(sigma_px,siz,pixels_per_mm,text)
    scale = siz./(1:siz)./pixels_per_mm;
    h=fspecial('gaussian', siz, sigma_px);
    hFT = abs(fft2(fftshift(h)));
    plot(scale,hFT(1,:)/max(hFT,[],'all'),'DisplayName',[text ' Gauss ' num2str(sigma_px/pixels_per_mm) 'mm'],'LineWidth',1.5)
    %imagesc(hFT/max(hFT,[],'all'))
end
