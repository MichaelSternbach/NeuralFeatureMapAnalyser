function runOPMDataTests(specimen_num_list,animal,bootstrapsamples,AnimalDataFolder,MainResultDataFolder)
    %runOPMDataTests(1:9,'dunnart',1000,'~/CIDBN/','~/Cloud/Cloud/PhD/MarsupialData/marsupial-data/ResultDatav1k/')
    
    %% set Figure file name
    FigureFileName = [MainResultDataFolder 'OPMDataTests_' animal];
    rm_cmd = ['rm -f ' FigureFileName '.ps'];
    disp(rm_cmd)
    system(rm_cmd)


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

    for specimen_num = specimen_num_list
        disp(['Processing specimen ', num2str(specimen_num)])

        %% Prepare a new figure for the current specimen
        h = 1300;
        f = figure('Visible', 'on', 'Position', [1, 1, h, 3*h]); % Adjust size for three vertically stacked subplots
        t = tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 3 rows, 1 column layout


        ResultDataFolder = [MainResultDataFolder lower(animal) '/' lower(animal) num2str(specimen_num) '/'];
        mkdir(ResultDataFolder)

        %% load data
        disp('load data')
        [data_info,~,data_obj,~,~] = getAnimalData(animal,specimen_num,AnimalDataFolder);
        data_obj.prepare_samples_array(bootstrapsamples)
        disp('-----------------------')

        %% get modularity distribution 
        DataFile = [ResultDataFolder 'ModularityDataMultiRand2.mat'];
        if ~isfile(DataFile)
            %$['~/Cloud/PhD/Writing/phd_thesis/OPM_Methods/Data/' animal '/' animal num2str(specimen_num) '/']
            [average_spacing_mm,~,~] =  getColumnsSpacing(data_obj,ResultDataFolder,smallest_w_mm,largest_w_mm,w_step_mm,false,FilterMap);
            peak_position_mm = average_spacing_mm;
            profile_range_mm_animal = sort([profile_range_mm average_spacing_mm]);
            [~,~,power_profiles,~,mean_abs_squared]=TestModularityOPM(data_obj,profile_range_mm_animal,profile_step_mm,Jackknife,FullCurve);
            [power_profiles_rand,mean_abs_squared_rand,~] =TestModularityOPM_MultiRand(data_obj,profile_range_mm_animal,profile_step_mm,bootstrapsamples);
            
            save(DataFile,'peak_position_mm',"power_profiles","power_profiles_rand",'mean_abs_squared','mean_abs_squared_rand','profile_range_mm_animal')
        else
            load(DataFile,'peak_position_mm',"power_profiles","power_profiles_rand",'mean_abs_squared','mean_abs_squared_rand','profile_range_mm_animal')
        end

%         %% plot power spectra
%         nexttile(t);
%         plotPowerspectra(power_profiles.BS,[0 1 1],'data')
%         plotPowerspectra(power_profiles_rand,[1 0 0],'rand')
%         hold on; plot([peak_position_mm peak_position_mm],[min(power_profiles.BS{1}.values_kspace) max(power_profiles.BS{1}.values_kspace)],'DisplayName','peak position')
%         legend()
%         xlabel('k [1/mm]')
%         ylabel('Power')
%         title(' Power Spectra')


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
        axis('square')
        title([animal num2str(specimen_num) ' Power Peak Distributions ' num2str(round(calcProbSmaller(peaks_test_rand,peaks_test(1)),3)) ' ' num2str(round(calcProbSmaller(peaks_test_rand,mean(peaks_test)),3)) ' ' num2str(round(calcProbSmaller(peaks_test_rand,median(peaks_test)),3))])   


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
        [x_high_prob_pw,y_high_prob_pw,Prob_high_prob_pw]=getPositionHighestProbPw(pinwheel_stats);
        selectivities_pw = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj);
        selectivities_pw_rand = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj_rand);
        
        %% plot histogram selectivity high prob. pinwheel
        nexttile(t);
        plotCPDFs(selectivities_pw,'Bootstrap Samples')
        hold on
        plotCPDFs(selectivities_pw_rand,'Randomized')
        hold on
        plot([selectivities_pw(1) selectivities_pw(1)],[0 1],'-red','DisplayName','mean map')
        hold on
        plot([mean(selectivities_pw) mean(selectivities_pw)],[0 1],'-','DisplayName','mean BS')
        hold on
        plot([median(selectivities_pw) median(selectivities_pw)],[0 1],'-','DisplayName','median BS')

        legend()
        xlabel('|z|^2')
        axis('square')
        title([ 'Pinwheel Highest Prob.' num2str(round(Prob_high_prob_pw,3)) ' Selectivties ' num2str(round(calcProbSmaller(selectivities_pw_rand,selectivities_pw(1)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand,mean(selectivities_pw)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand,median(selectivities_pw)),3))])
        %animal num2str(specimen_num)
        
        %% find sum selectivities pinwheels
        selectivities_pw_all = SumSelectivitiesPw(pinwheel_stats,data_obj);
        selectivities_pw_rand_all = SumSelectivitiesPw(pinwheel_stats,data_obj_rand);
        
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
        xlabel('\sum |z|^2')
        axis('square')
        title([animal num2str(specimen_num) 'Pinwheels Mean Selectivty ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,selectivities_pw_all(1)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,mean(selectivities_pw_all)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,median(selectivities_pw_all)),3))])


        %% Pinwheel Position Test
        linewidth = 0.2;
        ROI = data_obj.ROI;
        [YROI,XROI] = find(data_obj.ROI);
        [Xmin, Xmax] = findBorders(XROI);
        [Ymin, Ymax] = findBorders(YROI);

        nexttile(t);
        imagesc(zeros(size(ROI)))
        z = data_obj.filter_map(data_obj.read_map);
        plot_map(z,ROI,0,1)
        hold on;
        SizesCI = getConfidenceRegionPw(pinwheel_stats,data_info.field_size_pix,0.95,true,1000/data_info.pix_per_mm);
        
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
        plotCPDFs(sqrt(SizesCI)./data_info.pix_per_mm.*1000,'','-')
        title('Pinwheel CI Size CPDF')
        xlabel('$\sqrt{\text{PW CI size}} ≤ X [\mu m]$', 'Interpreter', 'latex')
        ylabel('% of pinwheels')
        %xlim([0,1])
        axis('square')

        %% Save the current specimen's plots to a PostScript file
        print(f, '-dpsc', '-fillpage', '-append', [FigureFileName, '.ps']);
        close(f); % Close the figure to save memory
    end
    %% Convert to pdf
    disp('Convert to pdf')
    cmd = ['ps2pdf ' FigureFileName '.ps ' FigureFileName, '.pdf'];
    disp(cmd)
    system(cmd)

end