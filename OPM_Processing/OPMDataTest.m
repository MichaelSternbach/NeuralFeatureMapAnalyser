function OPMDataTest(specimen_num,animal,bootstrapsamples,AnimalDataFolder,ResultDataFolder)
    % clear all
    % close all
    % disp('Start OPM analysis')
    % disp('-------------------------------------------')
    % %% add path to OPM processing functions 
    % addpath OPM_Processing/
    
%     
%     %% animal parameter
%     animal = 'dunnart';
%     %specimen_num = 1;
%     
%     bootstrapsamples = 100;
%     
%     %% data folder
%     AnimalDataFolder = '~/CIDBN/';#

    if ischar(specimen_num)
        specimen_num=str2num(specimen_num);
    end

    if ischar(bootstrapsamples)
        bootstrapsamples=str2num(bootstrapsamples);
    end
    
    disp([animal num2str(specimen_num)])
    %% parameter column spacing
    smallest_w_mm = 0.1;
    w_step_mm = 0.05;
    largest_w_mm = 1.5;
    FilterMap = true;
    
    %% parameter Testing
    
    Jackknife = false;
    
    profile_range_mm = [0.2 1.5];
    profile_step_mm = 0.01;
    profile_range_mm = profile_range_mm(1):profile_step_mm:profile_range_mm(2);
    FullCurve = true;
    
    % profile_range_mm = [0.01 1];
    % profile_step_mm = 0.01;
    % FullCurve = true;
    
    % for specimen_num = 1:9
%     ResultDataFolder = ['/home/michael/Cloud/PhD/MarsupialData/marsupial-data/Data/' lower(animal) '/' lower(animal) num2str(specimen_num) '/'];
    ResultDataFolder = [ResultDataFolder lower(animal) '/' lower(animal) num2str(specimen_num) '/'];
    mkdir(ResultDataFolder)

    %% load data
    disp('load data')
    [data_info,data_path,data_obj,data,BloodVesselImg] = getAnimalData(animal,specimen_num,AnimalDataFolder);
    data_obj.prepare_samples_array(bootstrapsamples)
    disp('-----------------------')
    %% get modularity distribution 

    DataFile = [ResultDataFolder 'ModularityDataMultiRand2.mat'];
    if ~isfile(DataFile)
        %$['~/Cloud/PhD/Writing/phd_thesis/OPM_Methods/Data/' animal '/' animal num2str(specimen_num) '/']
        [average_spacing_mm,local_spacing_mm,newROI] =  getColumnsSpacing(data_obj,ResultDataFolder,smallest_w_mm,largest_w_mm,w_step_mm,false,FilterMap);
        peak_position_mm = average_spacing_mm;
        profile_range_mm_animal = sort([profile_range_mm average_spacing_mm]);
        [~,~,power_profiles,~,mean_abs_squared]=TestModularityOPM(data_obj,profile_range_mm_animal,profile_step_mm,Jackknife,FullCurve);
        [power_profiles_rand,mean_abs_squared_rand,~] =TestModularityOPM_MultiRand(data_obj,profile_range_mm_animal,profile_step_mm,bootstrapsamples);
        
        save(DataFile,'peak_position_mm',"power_profiles","power_profiles_rand",'mean_abs_squared','mean_abs_squared_rand','profile_range_mm_animal')
    else
        load(DataFile,'peak_position_mm',"power_profiles","power_profiles_rand",'mean_abs_squared','mean_abs_squared_rand','profile_range_mm_animal')
    end
    
%     %% plot peak
%     f01=figure; plot(power_profiles.BS{1}.scale_mm,power_profiles.BS{1}.values)
%     hold on; plot([peak_position_mm pea   k_position_mm],[min(power_profiles.BS{1}.values) max(power_profiles.BS{1}.values)])
%     title([animal num2str(specimen_num)])
%     print(f01,'-depsc', [ResultDataFolder 'PeakPowerProfile.eps'])

    
    %% plot Power Hist
    fm1=figure;
    histogram(mean_abs_squared,'DisplayName','BS','FaceAlpha', 0.2)
    hold on;
    histogram(mean_abs_squared_rand,'DisplayName','Rand','FaceAlpha', 0.2)
    legend()
    xlabel('mean abs squared z')
    title([animal num2str(specimen_num) ' Power Hist'])
    print(fm1,'-depsc', [ResultDataFolder data_info.ID 'PowerHist.eps'])

% 
%     %% plot Modul profiles
%     f0=figure;
%     plotFullModularityTest(power_profiles.BS,[0 1 1],'data')
%     plotFullModularityTest(power_profiles_rand,[1 0 0],'rand')
%     if ~isnan(peak_position_mm)
%         hold on; plot([peak_position_mm peak_position_mm],[min(power_profiles.BS{1}.values) max(power_profiles.BS{1}.values)])
% %         hold on; plot([2*pi/peak_position_mm 2*pi/peak_position_mm],[min(power_profiles.BS{1}.values) max(power_profiles.BS{1}.values)])
%     end
%     %legend()
%     xlabel('r [mm]')
% %     xlabel('spacial frequnency [cycle/mm]')
%     ylabel('Power')
%     title([animal num2str(specimen_num) ' Power Profiles'])
% %     print(f0,'-depsc', [ResultDataFolder data_info.ID 'PowerProfilesMultiRand2.eps'])
%     print(f0, '-dpng', '-r300', [ResultDataFolder data_info.ID 'PowerProfilesMultiRand2.png']);
% 
%     
% 
%     
%     %% plot power spectra
%     f0=figure;
%     plotPowerspectra(power_profiles.BS,[0 1 1],'data')
%     plotPowerspectra(power_profiles_rand,[1 0 0],'rand')
%     if ~isnan(peak_position_mm)
% %         hold on; plot([peak_position_mm peak_position_mm],[min(power_profiles.BS{1}.values) max(power_profiles.BS{1}.values)])
%         hold on; plot([1/peak_position_mm 1/peak_position_mm],[min(power_profiles.BS{1}.values_kspace,[],'all') max(power_profiles.BS{1}.values_kspace,[],'all')])
%     end
%     xlabel('1/r [1/mm]')
% %     xlim([0.8 4])
%     ylabel('Power')
%     title([animal num2str(specimen_num) ' Power Spectra'])
% %     print(f0,'-depsc', [ResultDataFolder data_info.ID 'PowerSpectraMultiRand2.eps'])
%     print(f0, '-dpng', '-r300', [ResultDataFolder data_info.ID 'PowerSpectraMultiRand2.png'])
% 

    if ~isnan(peak_position_mm)
        %% check peak values
        I = find(power_profiles_rand{1}.scale_mm==peak_position_mm);
        peaks_test_rand = zeros([1 bootstrapsamples]);
        for ii = 1:bootstrapsamples
            peaks_test_rand(ii)=power_profiles_rand{ii}.values(I);
        end
    
        peaks_test = zeros([1 bootstrapsamples]);
        for ii = 1:bootstrapsamples
            peaks_test(ii)=power_profiles.BS{ii}.values(I);
        end
    
        f1=figure;
        n_bins = 15;
        max_peak = max([peaks_test peaks_test_rand],[],'all');
        min_peak = min([peaks_test peaks_test_rand],[],'all');
        bins = linspace(min_peak,max_peak,n_bins);
        %histogram(peaks_test,bins,'Normalization','pdf','DisplayName','Bootstrap Samples','FaceAlpha', 0.2)
        plotCPDFs(peaks_test,'Bootstrap Samples')
        hold on
%         histogram(peaks_test_rand,bins,'Normalization','pdf','DisplayName','Randomized','FaceAlpha', 0.2)
        plotCPDFs(peaks_test_rand,'Randomized')
        hold on
        plot([peaks_test(1) peaks_test(1)],[0 1],'-red','DisplayName','mean map')
        hold on
        plot([mean(peaks_test) mean(peaks_test)],[0 1],'-','DisplayName','mean BS')
        hold on
        plot([median(peaks_test) median(peaks_test)],[0 1],'-','DisplayName','median BS')
        legend()
        xlabel('peak hight')
        x=peaks_test;
        y=peaks_test_rand;
%         P_rand=sum(peaks_test_rand<=peaks_test(1),'all')/length(peaks_test_rand);
%         title([animal num2str(specimen_num) ' Power Peak Distributions ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5)) ' ' num2str(round(P_rand,5))])
        title([animal num2str(specimen_num) ' Power Peak Distributions ' num2str(round(calcProbSmaller(peaks_test_rand,peaks_test(1)),3)) ' ' num2str(round(calcProbSmaller(peaks_test_rand,mean(peaks_test)),3)) ' ' num2str(round(calcProbSmaller(peaks_test_rand,median(peaks_test)),3))])   
        print(f1,'-depsc', [ResultDataFolder data_info.ID 'ModularityDistributionMultiRand2.eps'])
    end


%     %% compare modularity distribution
%     n_bins = 15;
%     max_peak = max([peak_values peak_values_rand],[],'all');
%     min_peak = min([peak_values peak_values_rand],[],'all');
%     bins = linspace(min_peak,max_peak,n_bins);
%     
%     f1 = figure();
%     histogram(peak_values,bins,'Normalization','pdf','DisplayName','Bootstrap Samples')
%     hold on
%     histogram(peak_values_rand,bins,'Normalization','pdf','DisplayName','Randomized')
%     hold on
%     plot([peak_values(1) peak_values(1)],[0 0.12],'-red','DisplayName','mean map')
%     legend()
%     xlabel('peak hight')
%     title([animal num2str(specimen_num)])
%     print(f1,'-depsc', [ResultDataFolder 'ModularityDistribution' num2str(seed) '.eps'])

    %% get pinwheel stats
    tracker_obj = pinwheel_tracker;
    simple_track=true;
    
    DataFile = [ResultDataFolder 'PwStats.mat'];
    if ~isfile(DataFile)
        data_rand = randomizeData(data_obj.data);
        data_obj_rand =  data_handle_corrected(data_obj.info,data_rand,data_obj.ROI);
        data_obj.prepare_samples_array(bootstrapsamples)
        [pinwheel_stats,pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
    
    
        data_obj_rand.prepare_samples_array(bootstrapsamples)
        [pinwheel_stats_rand,pinwheel_spurious_rand] = get_pinwheel_stats(data_obj_rand,tracker_obj,simple_track);
    
    
        save(DataFile,'pinwheel_stats','pinwheel_spurious','pinwheel_stats_rand','pinwheel_spurious_rand','data_obj_rand')
    else
        load(DataFile,'pinwheel_stats','pinwheel_spurious','pinwheel_stats_rand','pinwheel_spurious_rand','data_obj_rand')
    end
    
    %% find selectivitieshighest prob. pinwheel
%     selectivities_pw = getSelectivitiesHighestProbPw(pinwheel_stats,data_obj);
%     selectivities_pw_rand = getSelectivitiesHighestProbPw(pinwheel_stats_rand,data_obj_rand);
    [x_high_prob_pw,y_high_prob_pw,Prob_high_prob_pw]=getPositionHighestProbPw(pinwheel_stats);
    selectivities_pw = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj);
    selectivities_pw_rand = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj_rand);
    
    %% plot histogram selectivity high prob. pinwheel
    n_bins = 20;
    max_peak = max([selectivities_pw;selectivities_pw_rand],[],'all');
    min_peak = min([selectivities_pw;selectivities_pw_rand],[],'all');
    bins = linspace(min_peak,max_peak,n_bins);
    
    f2 = figure;
%     histogram(selectivities_pw,bins,'Normalization','pdf','DisplayName','Bootstrap Samples','FaceAlpha', 0.2)
    plotCPDFs(selectivities_pw,'Bootstrap Samples')
    hold on
%     histogram(selectivities_pw_rand,bins,'Normalization','pdf','DisplayName','Randomized','FaceAlpha', 0.2)
    plotCPDFs(selectivities_pw_rand,'Randomized')
    hold on
    plot([selectivities_pw(1) selectivities_pw(1)],[0 1],'-red','DisplayName','mean map')
    hold on
    plot([mean(selectivities_pw) mean(selectivities_pw)],[0 1],'-','DisplayName','mean BS')
    hold on
    plot([median(selectivities_pw) median(selectivities_pw)],[0 1],'-','DisplayName','median BS')

    legend()
    xlabel('|z|^2')
    x=selectivities_pw;
    y=selectivities_pw_rand;
%     title([animal num2str(specimen_num) 'Pinwheel Highest Prob. Mean ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5))])
    title([animal num2str(specimen_num) 'Pinwheel Highest Prob.' num2str(round(Prob_high_prob_pw,3)) ' Selectivties ' num2str(round(calcProbSmaller(selectivities_pw_rand,selectivities_pw(1)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand,mean(selectivities_pw)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand,median(selectivities_pw)),3))])
    print(f2,'-depsc', [ResultDataFolder data_info.ID 'HighProbPwSelectivityDistribution2.eps'])
    
    
    %% find sum selectivities pinwheels
    selectivities_pw_all = SumSelectivitiesPw(pinwheel_stats,data_obj);
    selectivities_pw_rand_all = SumSelectivitiesPw(pinwheel_stats,data_obj_rand);
    
    %% plot histogram sum selectivity pinwheels
    
    
%     [chi2Stat, pValue, numBins] = compareDistributionsChiSquared(selectivities_pw_all, selectivities_pw_rand_all,5);
%     p = ranksum(selectivities_pw_all, selectivities_pw_rand_all);
%     p_width = RankSumWidth(selectivities_pw_all, selectivities_pw_rand_all);
    
    
    n_bins = 20;%numBins;
    max_peak = max([selectivities_pw_all;selectivities_pw_rand_all],[],'all');
    min_peak = min([selectivities_pw_all;selectivities_pw_rand_all],[],'all');
    bins = linspace(min_peak,max_peak,n_bins);
    
    f2 = figure;
%     histogram(selectivities_pw_all,bins,'Normalization','pdf','DisplayName','Bootstrap Samples','FaceAlpha', 0.2)
    plotCPDFs(selectivities_pw_all,'Bootstrap Samples')
    hold on
%     histogram(selectivities_pw_rand_all,bins,'Normalization','pdf','DisplayName','Randomized','FaceAlpha', 0.2)
    plotCPDFs(selectivities_pw_rand_all,'Randomized')
    hold on
    plot([selectivities_pw_all(1) selectivities_pw_all(1)],[0 1],'-red','DisplayName','mean map')
    hold on
    plot([mean(selectivities_pw_all) mean(selectivities_pw_all)],[0 1],'-','DisplayName','mean BS')
    hold on
    plot([median(selectivities_pw_all) median(selectivities_pw_all)],[0 1],'-','DisplayName','median BS')
    legend()
    xlabel('\sum |z|^2')
    x=selectivities_pw_all;
    y=selectivities_pw_rand_all;
%     title([animal num2str(specimen_num) 'Pinwheels Mean ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5))])
    title([animal num2str(specimen_num) 'Pinwheels Mean Selectivty ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,selectivities_pw_all(1)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,mean(selectivities_pw_all)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,median(selectivities_pw_all)),3))])
    
    print(f2,'-depsc', [ResultDataFolder data_info.ID 'SumPwSelectivityDistribution2.eps'])



    %% Pinwheel Position Test
    linewidth = 0.2;
    ROI = data_obj.ROI;
    [YROI,XROI] = find(data_obj.ROI);
    [Xmin, Xmax] = findBorders(XROI);
    [Ymin, Ymax] = findBorders(YROI);

    figure();
    imagesc(zeros(size(ROI)))
    z = data_obj.filter_map(data_obj.read_map);
    plot_map(z,ROI,0,1)
    hold on;
    SizesCI = getConfidenceRegionPw(pinwheel_stats,data_info.field_size_pix,0.95);
    
    contour(ROI,[1 1],'white','linewidth',linewidth)
%     m=100;
%     cm_viridis=viridis(m);
%     colormap( cm_viridis);
    plot(pinwheel_stats.x(:,1),pinwheel_stats.y(:,1),'wx')
    axis image
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
    title('95% CI Pinwheel Positions')
    %hold on; set(gca,'view',[rotate rotate])
    yticks([])
    xticks([])
    print('-depsc', [ResultDataFolder data_info.ID 'PwCI.eps'])

    
    
            
    %% plot CPDF pinwheel CI Size
    figure();
    PwCI = SizesCI/(data_info.pix_per_mm)^2;
    plotCPDFs(PwCI,'','-')
    title('Pinwheel CI Size CPDF')
    xlabel('PW CI size ≤ X [mm^2]')
    ylabel('% of pinwheels')
    %xlim([0,1])
    axis('square')
    print('-depsc', [ResultDataFolder data_info.ID 'PwCICPDF.eps'])
    
    close all
end






% ResultDataFolder = '/home/michael/Cloud/Cloud/PhD/MarsupialData/marsupial-data/ResultDatav1k/';
% AnimalDataFolder = '~/CIDBN/';
% for ii = 1:9
%     OPMDataTest(ii,'dunnart',1000,AnimalDataFolder,ResultDataFolder)
% end