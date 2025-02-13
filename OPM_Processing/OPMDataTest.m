function OPMDataTest(specimen_num_list,animal,bootstrapsamples,AnimalDataFolder,ResultDataFolder,GIF)

    if nargin < 6
        GIF = 4;
    end
    
    %% parameter Testing

    profile_range_mm = [0.1 2];
    profile_step_mm = 0.05;
    profile_range_mm = profile_range_mm(1):profile_step_mm:profile_range_mm(2);


    if ~isfolder(ResultDataFolder)
        mkdir(ResultDataFolder)
    end

    %% set filename
    FigureFileName = [ResultDataFolder animal 'TestOPM.ps'];

    %% remove old file
    if isfile(FigureFileName)
        cmd = ['rm ' FigureFileName];
        disp(cmd)
        disp(system(cmd))
    end


    for specimen_num = specimen_num_list

        
        disp([animal num2str(specimen_num)])
        ResultDataFolderAnimal = [ResultDataFolder lower(animal) '/' lower(animal) num2str(specimen_num) '/'];

        %% load data
        disp('load data')
        [data_info,~,data_obj,~,~] = getAnimalData(animal,specimen_num,AnimalDataFolder);
        data_obj.prepare_samples_array(bootstrapsamples)
        if GIF >0
            data_obj.activateGIF(true,GIF)
        end
        disp('-----------------------')

        [peak_position_mm,~,~] = getColumnsSpacing(data_obj,ResultDataFolderAnimal,profile_range_mm(1),profile_range_mm(2),profile_step_mm,false,false);
        
        %% get modularity distributions
        [peaks_test,peaks_test_rand,peak_position_mm,power_profiles,power_profiles_rand,~,~] = testModularityOPM(data_obj,ResultDataFolderAnimal,peak_position_mm,profile_range_mm,bootstrapsamples,false);

        %% get pw position data
        [pinwheel_stats,~,~,selectivities_pw_all,selectivities_pw_rand_all] = testPWsOPM(data_obj,[],bootstrapsamples,ResultDataFolderAnimal,false);
    

        %% set figure
        f = figure('Visible', 'on', 'Position', [1, 1, 800, 1600]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout
        title(t, [num2str(specimen_num) '. ' data_info.ID]);

        %% plot power spectra
        nexttile(t);
        plotPowerspectra(power_profiles.BS,[0 1 1],'data')
        plotPowerspectra(power_profiles_rand,[1 0 0],'rand')
        if ~isnan(peak_position_mm)
            hold on; 
            plot([1/peak_position_mm 1/peak_position_mm],[min(power_profiles.BS{1}.values_kspace) max(power_profiles.BS{1}.values_kspace)])
        end
        xlabel('k [1/mm]')
        ylabel('Power')
        title(' Power Spectra')

        %% plot Power Distribution
        nexttile(t);
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
        title(['Power Peak Distributions ' num2str(round(calcProbSmaller(peaks_test_rand,peaks_test(1)),3)) ' ' num2str(round(calcProbSmaller(peaks_test_rand,mean(peaks_test)),3)) ' ' num2str(round(calcProbSmaller(peaks_test_rand,median(peaks_test)),3))])   
        
        %% Pinwheel Position Test
        nexttile(t);
        linewidth = 0.2;
        ROI = data_obj.ROI;
        [YROI,XROI] = find(data_obj.ROI);
        [Xmin, Xmax] = findBorders(XROI);
        [Ymin, Ymax] = findBorders(YROI);

        imagesc(zeros(size(ROI)))
        z = data_obj.filter_map(data_obj.read_map);
        plot_map(z,ROI,0,1)
        hold on;
        SizesCI = getConfidenceRegionPw(pinwheel_stats,data_obj.info.field_size_pix,0.95,true,1000/data_obj.info.pix_per_mm);
        contour(ROI,[1 1],'white','linewidth',linewidth)
        plot(pinwheel_stats.x(:,1),pinwheel_stats.y(:,1),'wx')
        axis image
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('95% CI Pinwheel Positions')
        yticks([])
        xticks([])

        %% plot CPDF pinwheel CI Size
        nexttile(t);
        plotCPDFs(sqrt(SizesCI)./data_obj.info.pix_per_mm.*1000,'','-')
        title('Pinwheel CI Size CPDF')
        xlabel('sqrt(PW CI size) ≤ X [mu m]')
        ylabel('% of pinwheels')
        axis('square')

        %% plot histogram sum selectivity pinwheels
        nexttile(t);
        plotCPDFs(selectivities_pw_all,'Bootstrap Samples')
        hold on
        plotCPDFs(selectivities_pw_rand_all,'Randomized')
        hold on
        plot([selectivities_pw_all(1) selectivities_pw_all(1)],[0 1],'-red','DisplayName','mean map')
        hold on
        plot([mean(selectivities_pw_all) mean(selectivities_pw_all)],[0 1],'-','DisplayName','mean BS')
        hold on
        plot([median(selectivities_pw_all) median(selectivities_pw_all)],[0 1],'-','DisplayName','median BS')
        legend()
        xlabel('sum |z|^2')
        title(['Pinwheels Mean Selectivty ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,selectivities_pw_all(1)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,mean(selectivities_pw_all)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,median(selectivities_pw_all)),3))])
        


        %% save page and close figure
        print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
        close(f); % Close the figure to save memory
    end

end




