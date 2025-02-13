function [pinwheel_stats,selectivities_pw,selectivities_pw_rand,selectivities_pw_all,selectivities_pw_rand_all] = testPWsOPM(data_obj,pinwheel_stats,bootstrapsamples,ResultDataFolder,DoPlot)
    if nargin < 5
        DoPlot = true;
    end
    
    %% get pinwheel stats
    tracker_obj = pinwheel_tracker;
    simple_track=true;
    
    DataFile = [ResultDataFolder 'PwStats.mat'];
    if ~isfile(DataFile)

        if isempty(pinwheel_stats) || size(pinwheel_stats.x,2)~= bootstrapsamples
            disp('recalculate pinwheel stats')
            disp(size(pinwheel_stats.x,2))
            disp(bootstrapsamples)

            data_obj.prepare_samples_array(bootstrapsamples)
            [pinwheel_stats,~] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
            
        end
        
        %% prepare randomized data
        data_rand = randomizeData(data_obj.data);
        data_obj_rand =  data_handle_corrected(data_obj.info,data_rand,data_obj.ROI);
        if data_obj.GIF_apply
            data_obj_rand.activateGIF(true,data_obj.SN_TH)
        end
        if data_obj.lsm_applied
            data_obj_rand.apply_LSM(true)
        end
        data_obj_rand.prepare_samples_array(bootstrapsamples)

    
        %% find selectivitieshighest prob. pinwheel
        [x_high_prob_pw,y_high_prob_pw,Prob_high_prob_pw]=getPositionHighestProbPw(pinwheel_stats);
        selectivities_pw = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj);
        selectivities_pw_rand = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj_rand);

        %% find sum selectivities pinwheels
        selectivities_pw_all = SumSelectivitiesPw(pinwheel_stats,data_obj);
        selectivities_pw_rand_all = SumSelectivitiesPw(pinwheel_stats,data_obj_rand);
    
        save(DataFile,'pinwheel_stats','data_obj_rand', ...
            'Prob_high_prob_pw','selectivities_pw','selectivities_pw_rand',"selectivities_pw_all",'selectivities_pw_rand_all')
    else
        load(DataFile,'pinwheel_stats', ...
            'Prob_high_prob_pw','selectivities_pw','selectivities_pw_rand',"selectivities_pw_all",'selectivities_pw_rand_all')
    end
    
    if DoPlot
        %% plot distribution selectivity high prob. pinwheel
        
        f2 = figure;
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
        title([data_obj.info.animal ' ' data_obj.info.ID 'Pinwheel Highest Prob.' num2str(round(Prob_high_prob_pw,3)) ' Selectivties ' num2str(round(calcProbSmaller(selectivities_pw_rand,selectivities_pw(1)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand,mean(selectivities_pw)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand,median(selectivities_pw)),3))])
        print(f2,'-depsc', [ResultDataFolder data_obj.info.ID 'HighProbPwSelectivityDistribution2.eps'])
        
        

        
        %% plot histogram sum selectivity pinwheels
        f2 = figure;
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
        title([data_obj.info.animal ' ' data_obj.info.ID 'Pinwheels Mean Selectivty ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,selectivities_pw_all(1)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,mean(selectivities_pw_all)),3)) ' ' num2str(round(calcProbSmaller(selectivities_pw_rand_all,median(selectivities_pw_all)),3))])
        
        print(f2,'-depsc', [ResultDataFolder data_obj.info.ID 'SumPwSelectivityDistribution2.eps'])



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
        SizesCI = getConfidenceRegionPw(pinwheel_stats,data_obj.info.field_size_pix,0.95,true,1000/data_obj.info.pix_per_mm);
        contour(ROI,[1 1],'white','linewidth',linewidth)
        plot(pinwheel_stats.x(:,1),pinwheel_stats.y(:,1),'wx')
        axis image
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        title('95% CI Pinwheel Positions')
        yticks([])
        xticks([])
        print('-depsc', [ResultDataFolder data_obj.info.ID 'PwCI.eps'])

        
        
                
        %% plot CPDF pinwheel CI Size
        figure();
        plotCPDFs(sqrt(SizesCI)./data_obj.info.pix_per_mm.*1000,'','-')
        title('Pinwheel CI Size CPDF')
        xlabel('sqrt(PW CI size) ≤ X [mu m]')
        ylabel('% of pinwheels')
        axis('square')
        print('-depsc', [ResultDataFolder data_obj.info.ID 'PwCICPDF.eps'])
        
        close all
    end
end