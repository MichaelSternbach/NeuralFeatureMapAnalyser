function testModularityOPM(data_obj,ResultDataFolder,peak_position_mm,profile_range_mm,bootstrapsamples)
    
    Jackknife = false;
    FullCurve = true;
    
    DataFile = [ResultDataFolder 'TestModularityOPM.mat'];
    if ~isfile(DataFile)
        profile_range_mm_animal = sort(unique([profile_range_mm peak_position_mm]));
        
        data_obj.prepare_samples_array(bootstrapsamples)
        [~,~,power_profiles,~,mean_abs_squared] = TestModularityOPM(data_obj,profile_range_mm_animal,Jackknife,FullCurve);
        
        [power_profiles_rand,mean_abs_squared_rand,~] =TestModularityOPM_MultiRand(data_obj,profile_range_mm_animal,bootstrapsamples);
        
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
    title([' Power Hist'])
    print(fm1,'-depsc', [ResultDataFolder data_obj.info.ID 'PowerHist.eps'])


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
    title([ ' Power Profiles'])
    print(f0,'-depsc', [ResultDataFolder data_obj.info.ID 'PowerProfilesMultiRand.eps'])
    

    
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
    title([ ' Power Spectra'])
    print(f0,'-depsc', [ResultDataFolder data_obj.info.ID 'PowerSpectraMultiRand.eps'])

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
        title([ ' Power Peak Distributions ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5))])
        print(f1,'-depsc', [ResultDataFolder data_obj.info.ID 'ModularityDistributionMultiRand.eps'])
    end
end