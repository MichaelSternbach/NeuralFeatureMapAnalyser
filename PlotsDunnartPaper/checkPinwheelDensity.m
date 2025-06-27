function checkPinwheelDensity(experiment_num_list,animal,AnimalDataFolder,DataFolderMain,...
    GIF_SNTH,FigureFolder)
    %checkPinwheelDensity(1:6,'wallaby','~/CIDBN/','/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/DataHPC_GIF/')
    close all
    %% parameter
    experiment_num_list = checkFormatNum(experiment_num_list);
    if nargin <5
        GIF_SNTH = 4;
    end
    if nargin<6
        FigureFolder = DataFolderMain;
    end
    width_scale_pix = 15;
    linewidth = 0.3;
    smallest_w_mm = 0.4;
    largest_w_mm=2;
    w_step_mm=0.1;
%     
%     %% set filter settings
%     lowpass_mm = 0.25;
%     highpass_mm = 0.56;
    
    %% figure folder and files
    if ~isfolder(FigureFolder)
        mkdir(FigureFolder) 
    end

    FigureFile1 = [FigureFolder animal '_PinwheelsGIF' num2str(GIF_SNTH)];    
    FigureFile2 = [FigureFolder animal '_ColumnSpacingGIF' num2str(GIF_SNTH)];
    
    %% make Figures
    N_Animals = length(experiment_num_list);

    f1 = figure(1);
    t1 = tiledlayout(ceil(N_Animals/3),3);
    %title(t1,'PwDesnity Kernel/PD Simple/#Pw/#Columns')

    f2 = figure(2);
    t2 = tiledlayout(ceil(N_Animals/3),3);
    %title(t2,'PwDesnity Kernel/PD Simple/#Pw/#Columns')

    

    T = table();
    for experiment_num = experiment_num_list
        
        %% data folder
        DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
        
        %% animal data
        [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,AnimalDataFolder);
        if GIF_SNTH>0
            data_obj.activateGIF(true,GIF_SNTH)
        end

        %% rectangle ROI
        rectangle = data_info.rectangle;
        rectangleROI = zeros(size(data_obj.ROI));
        rectangleROI(rectangle(2):rectangle(4),rectangle(1):rectangle(3))=1;
        data_obj.set_ROI(rectangleROI);

        %% make map borders ROI
        [YROI,XROI] = find(data_obj.ROI);
        [Xmin, Xmax] = findBorders(XROI);
        [Ymin, Ymax] = findBorders(YROI);
        
        %% get map
        z = data_obj.filter_map(data_obj.read_map());
        
%         im = imag(z);
%         re = real(z);
% 
%         im = im/std(im,0,'all');
%         re = re/std(re,0,'all');
% 
%         z = re + 1i*im;


        %% find pinwheels Jason
        [pinw, windno] = locate_pinwheels(z.*data_obj.ROI);

        %% get column spacing Jason
        op_wl = mean([fourier_wavelength(real(z.*data_obj.ROI)), fourier_wavelength(imag(z.*data_obj.ROI))]);

        %% calc pinwheel density Jason
        npinw = size(pinw,1);

        % Calculate per pixel
        density_per_pixel = npinw/sum(data_obj.ROI(:));
        
        % Calculate per wavelength
        density_per_wl = density_per_pixel*op_wl^2;

        op_wl_mm = op_wl/data_info.pix_per_mm;
        
        %% calculate column spacing
        [average_spacing_mmManuel,local_spacing_mmManuel,~,WavletCoefficient] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
        [average_spacing_mm,local_spacing_mm] = get_column_spacing(z,data_obj.ROI,data_info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm,0);
        PwInfo.average_spacing_mm = average_spacing_mmManuel;
        PwInfo.NumHypercolumns = sum(data_obj.ROI,'all')/(data_obj.info.pix_per_mm*PwInfo.average_spacing_mm)^2;
        
        %% find pinwheels
        [PwInfo.NumberPw,PwInfo.aniso,PwInfo.x_angle,PwInfo.PWxList,PwInfo.PWyList,PwInfo.signList, PwInfo.contours] = find_pinwheels(z,0,data_obj.ROI);
        

        %% get Local PwDensity Fixed Filter
        sigma = 0.1;
        PwInfo.LocalPwDensityFixedFilter=getLocalPwDensityFixedFilter(data_obj,PwInfo,local_spacing_mmManuel,sigma);
        PwInfo.WeightedPwDensityFixedFilter = mean(PwInfo.LocalPwDensityFixedFilter(data_obj.ROI));

        %% get simple pw density 
        PwInfo.MeanPwDensity = PwInfo.NumberPw/PwInfo.NumHypercolumns;


        f1 = figure(1);
        %% plot map
        nexttile(t1)
        plot_map(z,data_obj.ROI,0,1)
        
        %% plot column spacing
        scale = (data_info.pix_per_mm*average_spacing_mmManuel).^-1;
        spacing_pix = 1/scale;
        hold on
        plot([width_scale_pix,width_scale_pix+spacing_pix]+Xmin,[width_scale_pix, width_scale_pix]+Ymin,'-white')
        hold on
        text(Xmin+width_scale_pix+spacing_pix+2,width_scale_pix+Ymin,'Λ','Color','white')

        hold on
        plot([width_scale_pix,width_scale_pix+op_wl]+Xmin,[2*width_scale_pix, 2*width_scale_pix]+Ymin,'--white')
        
        %% plot pinwheels
        hold on; plot(PwInfo.PWxList,PwInfo.PWyList,'o','Color','white','linewidth',linewidth)

        %% plot pinwheels Jason
        hold on; plot(pinw(:,1),pinw(:,2),'x','Color','white','linewidth',linewidth)
        

        %% add ROI limits
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])

        %% plot pw densities and numbers
%         title([num2str(PwInfo.WeightedPwDensityFixedFilter) '/' num2str(PwInfo.MeanPwDensity) '/' num2str(PwInfo.NumberPw) '/' num2str(PwInfo.NumHypercolumns)])
        title([ num2str(PwInfo.MeanPwDensity) '/' num2str(PwInfo.NumberPw) '/' num2str(average_spacing_mmManuel) 'mm ||' num2str(density_per_wl) '/' num2str(npinw) '/' num2str(op_wl_mm) 'mm'])
        

        %% calculate power profile 
        power_spectrum_range_mm = smallest_w_mm:0.05:largest_w_mm;
        power_profile_filtered = define_filter_settings(data_info,data_obj.ROI,z,power_spectrum_range_mm);
        power_profile = define_filter_settings(data_info,data_obj.ROI,data_obj.read_map(),power_spectrum_range_mm);
        
        
        %% column spacing
        f2 = figure(2);
        nexttile(t2)
        plot(WavletCoefficient.X, WavletCoefficient.Y_mean/max(WavletCoefficient.Y_mean),'DisplayName','mean wavelet');
        hold on
        plot(WavletCoefficient.XI, WavletCoefficient.YI_mean/max(WavletCoefficient.Y_mean),'DisplayName','mean wavelet fit');
        hold on
        plot(power_profile_filtered.scale_mm*1000,power_profile_filtered.values/max(power_profile_filtered.values),'DisplayName','power profile filtered');
        hold on
        plot(power_profile.scale_mm*1000,power_profile.values/max(power_profile.values),'DisplayName','power profile');
        hold on
        plot([average_spacing_mm*1000 average_spacing_mm*1000],[0 1],'-','DisplayName','Wavelet symmetric');
        hold on
        plot([average_spacing_mmManuel*1000 average_spacing_mmManuel*1000],[0 1],'DisplayName','Wavelet gabor');
        hold on
        plot([op_wl_mm*1000 op_wl_mm*1000],[0 1],'DisplayName','Fourier');
        title([num2str(experiment_num) '. ' data_info.ID])
        xlabel('Scale [mu m]');
        legend();
       % ylabel('Wavelet Coeff');

        %% save results to table
        T = [T; getRowPwTable(data_info,PwInfo)];
    end
    
    %% save Plots
    print(f1, '-dpsc','-fillpage', [FigureFile1 '.ps'])%'-append',
    print(f2, '-dpsc','-fillpage', [FigureFile2 '.ps'])%'-append',

    %% save figures
    saveas(f1, [FigureFile1 '.fig'])
    saveas(f2, [FigureFile2 '.fig'])

    %% save Table
%     writetable(T, [FigureFile1 '.txt'], '\t');
    
    
    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return 
end


function Row = getRowPwTable(data_info,PwInfo)
   Name = string(data_info.ID);
%    PwD = num2str(PwInfo.WeightedPwDensityFixedFilter);
%    num_HC = num2str(PwInfo.NumHypercolumns);
%    %Age_months = num2str(data_info.age);
%    Weight_g = num2str(data_info.weight_in_grams);
%    PWL = num2str(PwInfo.NumberPw);
%    WL_mm = num2str(PwInfo.average_spacing_mm);
   PwD = PwInfo.WeightedPwDensityFixedFilter;
   num_HC = PwInfo.NumHypercolumns;
   %Age_months = data_info.age);
   Weight_g = data_info.weight_in_grams;
   PWL = PwInfo.NumberPw;
   WL_mm = PwInfo.average_spacing_mm;

%    Row = table(Name, PwD, num_HC, Age_months,Weight_g,PWL,WL_mm); 
    Row = table(Name, PwD, num_HC,Weight_g,PWL,WL_mm); 
end
    


function LocalPwDensityFixedFilter = getLocalPwDensityFixedFilter(data_obj,PwInfo,local_spacing_mm,sigma)
    disp('calc LocalPwDensity with FixedFilter')
    average_spacing_mm = mean(local_spacing_mm(data_obj.ROI));
    local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PwInfo.PWxList, PwInfo.PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
    local_pw_dens = local_pw_dens./sum(local_pw_dens(data_obj.ROI)).*PwInfo.NumberPw;
    LocalPwDensityFixedFilter = local_pw_dens.*(local_spacing_mm*data_obj.info.pix_per_mm).^2;
end

