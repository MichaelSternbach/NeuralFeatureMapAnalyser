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
    print(fm1,'-depsc', [ResultDataFolder 'PowerHist.eps'])


    %% plot Modul profiles
    f0=figure;
    plotFullModularityTest(power_profiles.BS,[0 1 1],'data')
    plotFullModularityTest(power_profiles_rand,[1 0 0],'rand')
    if ~isnan(peak_position_mm)
        hold on; plot([peak_position_mm peak_position_mm],[min(power_profiles.BS{1}.values) max(power_profiles.BS{1}.values)])
    end
    %legend()
    xlabel('r [mm]')
    ylabel('Power')
    title([animal num2str(specimen_num) ' Power Profiles'])
    print(f0,'-depsc', [ResultDataFolder 'PowerProfilesMultiRand2.eps'])
    

    
    %% plot power spectra
    f0=figure;
    plotPowerspectra(power_profiles.BS,[0 1 1],'data')
    plotPowerspectra(power_profiles_rand,[1 0 0],'rand')
%     if ~isnan(peak_position_mm)
%         hold on; plot([peak_position_mm peak_position_mm],[min(power_profiles.BS{1}.values_kspace) max(power_profiles.BS{1}.values_kspace)])
%     end
    %legend()
    xlabel('k [1/mm]')
    ylabel('Power')
    title([animal num2str(specimen_num) ' Power Spectra'])
    print(f0,'-depsc', [ResultDataFolder 'PowerSpectraMultiRand2.eps'])

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
        histogram(peaks_test,bins,'Normalization','pdf','DisplayName','Bootstrap Samples','FaceAlpha', 0.2)
        hold on
        histogram(peaks_test_rand,bins,'Normalization','pdf','DisplayName','Randomized','FaceAlpha', 0.2)
        hold on
        plot([peaks_test(1) peaks_test(1)],[0 0.12],'-red','DisplayName','mean map')
        legend()
        xlabel('peak hight')
        x=peaks_test;
        y=peaks_test_rand;
        title([animal num2str(specimen_num) ' Power Peak Distributions ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5))])
        print(f1,'-depsc', [ResultDataFolder 'ModularityDistributionMultiRand2.eps'])
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
    [x_high_prob_pw,y_high_prob_pw]=getPositionHighestProbPw(pinwheel_stats);
    selectivities_pw = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj);
    selectivities_pw_rand = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj_rand);
    
    %% plot histogram selectivity high prob. pinwheel
    n_bins = 20;
    max_peak = max([selectivities_pw;selectivities_pw_rand],[],'all');
    min_peak = min([selectivities_pw;selectivities_pw_rand],[],'all');
    bins = linspace(min_peak,max_peak,n_bins);
    
    f2 = figure;
    histogram(selectivities_pw,bins,'Normalization','pdf','DisplayName','Bootstrap Samples','FaceAlpha', 0.2)
    hold on
    histogram(selectivities_pw_rand,bins,'Normalization','pdf','DisplayName','Randomized','FaceAlpha', 0.2)
    
    legend()
    xlabel('|z|^2')
    x=selectivities_pw;
    y=selectivities_pw_rand;
    title([animal num2str(specimen_num) 'Pinwheel Highest Prob. Mean ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5))])
    print(f2,'-depsc', [ResultDataFolder 'HighProbPwSelectivityDistribution2.eps'])
    
    
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
    histogram(selectivities_pw_all,bins,'Normalization','pdf','DisplayName','Bootstrap Samples','FaceAlpha', 0.2)
    hold on
    histogram(selectivities_pw_rand_all,bins,'Normalization','pdf','DisplayName','Randomized','FaceAlpha', 0.2)
    
    legend()
    xlabel('\sum |z|^2')
    x=selectivities_pw_all;
    y=selectivities_pw_rand_all;
    title([animal num2str(specimen_num) 'Pinwheels Mean ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5))])
    print(f2,'-depsc', [ResultDataFolder 'SumPwSelectivityDistribution2.eps'])
    
    close all
end