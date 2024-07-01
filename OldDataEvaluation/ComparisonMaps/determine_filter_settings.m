for ferret = 1:1
    llp_cutoffs = linspace(0.2,1,0.2);
    
    
    clf%figure('color','w')
    
    disp(ferret)
    
    [data_info,data_path] = info_handle('ferret',ferret);
    trials_to_use = find(data_info.protocol.blocks>0);
    
    max_plot = 0;
    prc_sum = zeros(size(llp_cutoffs));
    for trial = trials_to_use
        
        load([data_path,'/Analyzed_new/characterization/trial_',num2str(trial),'_design.mat'],'design','spectrum_profile','spectrum_scale_mm','llp_cutoffs')
         max_plot = max([max_plot max(spectrum_profile)]);
        % radial profile
        subplot(1,2,1)
        hold on
         plot(spectrum_scale_mm,spectrum_profile,'k')
         hold on
         ind = find(spectrum_scale_mm<=design.average_w_microns/1000,1);
         plot(spectrum_scale_mm(ind),spectrum_profile(ind),'ok','markerfacecolor','r')
         
          % plateau
        pixels_per_mm = design.measure;
        local_w = design.local_w(design.new_roi)/pixels_per_mm;
        tmp = design.cutoff_lowpass(:,:,1);
        lp_low = tmp(design.new_roi);
        tmp = design.cutoff_lowpass(:,:,2);
        lp_high = tmp(design.new_roi);
                
        prc = zeros(size(llp_cutoffs));
        for ii=1:length(llp_cutoffs)
            test = llp_cutoffs(ii)>=lp_low & llp_cutoffs(ii)<=lp_high;
            prc(ii) = 100*mean(test);
        end
        
        subplot(1,2,2)
        hold on
        plot(llp_cutoffs,prc,'k')
        
        prc_sum = prc_sum + prc/length(trials_to_use);
    end
    
    % get max
    [~,ind] = max(prc_sum);
    lp_to_use = llp_cutoffs(ind);
    hp_to_use = 1.5;
    
    % labels
    subplot(1,2,1)
    hold on
    line([lp_to_use lp_to_use],[0 max_plot*1.1],'color','g')
    line([hp_to_use hp_to_use],[0 max_plot*1.1],'color','g')
    xlabel('Scale in mm')
    ylabel('Power')
    xlim([0 max(spectrum_scale_mm)])
    title('Spectral profile')
    
    subplot(1,2,2)
    plot(llp_cutoffs,prc_sum,'r')
    line([lp_to_use lp_to_use],[0 100],'color','g')
    xlabel('LP cutoff')
    ylabel('% in plateau')
    ylim([0 100])
    title(['LP = ',num2str(lp_to_use)])
    
    export_fig([data_path,'Processed_GIF/cutoff_vs_plateau.jpg'],'-djpg','-m2','-nocrop');
end

%%

%%
for ferret = 1:31
    [data_info,data_path] = info_handle('ferret',ferret);
    copyfile([data_path,'Processed_GIF/cutoff_vs_plateau.jpg'],['/pairing/PairingData/Ferret_Whitney/Analysis/filter/',data_info.ID,'.jpg']);
end

%%
for ferret = 1:31
    [data_info,data_path] = info_handle('ferret',ferret);
    copyfile([data_path,'/Processed_GIF/layouts_filtered.jpg'],['/pairing/PairingData/Ferret_Whitney/Analysis/filter/',data_info.ID,'_layout.jpg']);
end
