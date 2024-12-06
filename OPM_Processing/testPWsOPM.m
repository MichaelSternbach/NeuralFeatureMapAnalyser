function testPWsOPM(data_obj,pinwheel_stats,pinwheel_spurious,bootstrapsamples,ResultDataFolder)
    %% get pinwheel stats
    tracker_obj = pinwheel_tracker;
    simple_track=true;
    
    DataFile = [ResultDataFolder 'PwStats.mat'];
    if ~isfile(DataFile)

        if size(pinwheel_stats.x,2) ~= bootstrapsamples
            disp('recalculate pinwheel stats')
            disp(size(pinwheel_stats.x,2))
            disp(bootstrapsamples)

            data_obj.prepare_samples_array(bootstrapsamples)
            [pinwheel_stats,pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
            
        end

        data_rand = randomizeData(data_obj.data);
        data_obj_rand =  data_handle_corrected(data_obj.info,data_rand,data_obj.ROI);    
        data_obj_rand.prepare_samples_array(bootstrapsamples)
        [pinwheel_stats_rand,pinwheel_spurious_rand] = get_pinwheel_stats(data_obj_rand,tracker_obj,simple_track);

%         data_obj.prepare_samples_array(bootstrapsamples)
%         [pinwheel_stats,pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
    
        %% find selectivitieshighest prob. pinwheel
%         selectivities_pw = getSelectivitiesHighestProbPw(pinwheel_stats,data_obj);
%         selectivities_pw_rand = getSelectivitiesHighestProbPw(pinwheel_stats_rand,data_obj_rand);
        [x_high_prob_pw,y_high_prob_pw]=getPositionHighestProbPw(pinwheel_stats);
        selectivities_pw = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj);
        selectivities_pw_rand = getSelectivitiesPw(x_high_prob_pw,y_high_prob_pw,data_obj_rand);

        %% find sum selectivities pinwheels
        selectivities_pw_all = SumSelectivitiesPw(pinwheel_stats,data_obj);
        selectivities_pw_rand_all = SumSelectivitiesPw(pinwheel_stats,data_obj_rand);
    
        save(DataFile,'pinwheel_stats','pinwheel_spurious','pinwheel_stats_rand','pinwheel_spurious_rand','data_obj_rand', ...
            'selectivities_pw','selectivities_pw_rand',"selectivities_pw_all",'selectivities_pw_rand_all')
    else
        load(DataFile,'pinwheel_stats','pinwheel_spurious','pinwheel_stats_rand','pinwheel_spurious_rand','data_obj_rand', ...
            'selectivities_pw','selectivities_pw_rand',"selectivities_pw_all",'selectivities_pw_rand_all')
    end
    

    
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
    xlabel('|z_pw|^2')
    x=selectivities_pw;
    y=selectivities_pw_rand;
    title([ 'Pinwheel Highest Prob. Mean ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5))])
    print(f2,'-depsc', [ResultDataFolder data_obj.info.ID 'HighProbPwSelectivityDistribution.eps'])
    
    

    
    %% plot histogram sum selectivity pinwheels 
    
    n_bins = 20;%numBins;
    max_peak = max([selectivities_pw_all;selectivities_pw_rand_all],[],'all');
    min_peak = min([selectivities_pw_all;selectivities_pw_rand_all],[],'all');
    bins = linspace(min_peak,max_peak,n_bins);
    
    f2 = figure;
    histogram(selectivities_pw_all,bins,'Normalization','pdf','DisplayName','Bootstrap Samples','FaceAlpha', 0.2)
    hold on
    histogram(selectivities_pw_rand_all,bins,'Normalization','pdf','DisplayName','Randomized','FaceAlpha', 0.2)
    
    legend()
    xlabel('sum |z_pw|^2')
    x=selectivities_pw_all;
    y=selectivities_pw_rand_all;
    title([ 'Pinwheels Mean ' num2str(round(ranksum(x,y),5)) ' ' num2str(round(RankSumWidth(x,y),5))])
    print(f2,'-depsc', [ResultDataFolder data_obj.info.ID 'SumPwSelectivityDistribution.eps'])

    close all
end