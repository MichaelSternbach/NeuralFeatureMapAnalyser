function [peaks_test,peaks_test_rand,peak_position_mm,power_profiles,power_profiles_rand,mean_abs_squared,mean_abs_squared_rand] = testModularityOPM(data_obj,ResultDataFolder,peak_position_mm,profile_range_mm,bootstrapsamples,DoPlot)

    if nargin < 6
        DoPlot = true;
    end
    
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
        load(DataFile,'peak_position_mm',"power_profiles","power_profiles_rand",'mean_abs_squared','mean_abs_squared_rand')
    end

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
    end

    if DoPlot
        %% plot Power Hist
        fm1=figure;
        histogram(mean_abs_squared,'DisplayName','BS','FaceAlpha', 0.2)
        hold on;
        histogram(mean_abs_squared_rand,'DisplayName','Rand','FaceAlpha', 0.2)
        legend()
        xlabel('mean abs squared z')
        title([data_obj.info.animal ' ' data_obj.info.ID ' Power Hist'])
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
        title(' Power Profiles')
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
        title(' Power Spectra')
        print(f0,'-depsc', [ResultDataFolder data_obj.info.ID 'PowerSpectraMultiRand.eps'])

        if ~isnan(peak_position_mm)
        
            f1=figure;
            plotCPDFs(peaks_test,'Bootstrap Samples')
            hold on
            plotCPDFs(peaks_test_rand,'Randomized')
            hold on
            plot([peaks_test(1) peaks_test(1)],[0 1],'-red','DisplayName','mean map')
            hold on
            plot([mean(peaks_test) mean(peaks_test)],[0 1],'-','DisplayName','mean BS')
            hold on
            plot([median(peaks_test) median(peaks_test)],[0 1],'-','DisplayName','median BS')
            legend()
            xlabel('peak hight')
            title([data_obj.info.animal ' ' data_obj.info.ID ' Power Peak Distributions ' num2str(round(calcProbSmaller(peaks_test_rand,peaks_test(1)),3)) ' ' num2str(round(calcProbSmaller(peaks_test_rand,mean(peaks_test)),3)) ' ' num2str(round(calcProbSmaller(peaks_test_rand,median(peaks_test)),3))])   
            print(f1,'-depsc', [ResultDataFolder data_obj.info.ID 'ModularityDistributionMultiRand2.eps'])
        end
    end
end