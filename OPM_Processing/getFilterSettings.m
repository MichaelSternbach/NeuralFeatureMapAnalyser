function data_info = getFilterSettings(data_obj,data_info,average_spacing_mm,folder,resetFilter, lowpass_cutoffs,profile_range_mm,profile_step_mm)
    if nargin <5
        resetFilter = false;
    end
    if nargin <6
        lowpass_cutoffs = linspace(0.2*average_spacing_mm,average_spacing_mm,50);
    end
    
    if nargin <7
        profile_range_mm = [0.01 2];
    end
    if nargin <8 
        profile_step_mm = 0.01;
    end
    
    DataAndFilterFile = [folder data_info.ID '.mat'];
    
    %% input filter parameter
    if resetFilter || ~isfield(data_info,'settings') || ~isfield(data_info.settings,'lowpass_mm') || ~isfield(data_info.settings,'highpass_mm')
        power_profile = define_filter_settings(data_info,data_obj.ROI,data_obj.data,profile_range_mm,profile_step_mm);

        figure
        plot(power_profile.scale_mm,power_profile.values);
        xlabel('Scale in mm')
        ylabel('Power')
        xlim([0 max(profile_range_mm)])
        set(gca,'fontsize',15)

        prompt = "What is the lowpass filter value in mm? ";
        lowpass_mm = input(prompt);
        data_obj.set_filter_parameters('lowpass',lowpass_mm)
        data_info.settings.lowpass_mm = lowpass_mm;

        prompt = "What is the highpass filter value in mm? ";
        highpass_mm = input(prompt);
        data_obj.set_filter_parameters('highpass',highpass_mm)
        data_info.settings.highpass_mm = highpass_mm;

        close all


        %% determine lowpass filter from pinwheel plateau

        lowpass_cutoffs = lowpass_cutoffs(lowpass_cutoffs<(highpass_mm*.9));
        filtersPwNumber = find_lowpassPwNumber(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
        
        f=figure;
        plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2),'DisplayName','pinwheel data')
        xlabel('lowpassfilter cutoff [mm]')
        ylabel('pinwheel number')

        prompt = "What is the lowpass filter value in mm based on the pinwheel density plateau? ";
        lowpass_mm = input(prompt);
        data_obj.set_filter_parameters('lowpass',lowpass_mm)
        data_info.settings.lowpass_mm = lowpass_mm;
        close(f)

        %% Last Check Filter parameter

        figure

        plot(power_profile.scale_mm,power_profile.values,'DisplayName','power profile unfiltered map');
        hold on
        plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','lowpass cutoff')
        hold on
        plot([data_obj.filter_parameters.highpass data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','highpass cutoff')

        xlabel('Scale in mm')
        ylabel('Power')
        xlim([0.01 2])
        ylim([min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
        set(gca,'fontsize',15)
        legend()
        
        
        figure 
        plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2),'DisplayName','pinwheel data')
        hold on
        plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all')],'DisplayName','lowpass cutoff')

        xlabel('Scale in mm')
        ylabel('# pinwheels')
        set(gca,'fontsize',15)
        legend()


        disp(['lowpassfilter cutoff [mm]: ' num2str(data_obj.filter_parameters.lowpass)])
        disp(['highpassfilter cutoff [mm]: ' num2str(data_obj.filter_parameters.highpass)])

        prompt = "Are the filter parameter right?  ";
        Answer = char(lower(input(prompt)));

        if Answer(1) ~= 'y' 
            prompt = "What is the lowpass filter value in mm? ";
            lowpass_mm = input(prompt);
            data_obj.set_filter_parameters('lowpass',lowpass_mm)
            data_info.settings.lowpass_mm = lowpass_mm;

            prompt = "What is the highpass filter value in mm? ";
            highpass_mm = input(prompt);
            data_obj.set_filter_parameters('highpass',highpass_mm)
            data_info.settings.highpass_mm = highpass_mm;
        end

        close all

        %% save (filter)parameter and data
        save_info(data_obj.info.data_path,data_info)
        save(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumber','power_profile')
    else
        if isfile(DataAndFilterFile)
            load(DataAndFilterFile,'filtersPwNumber','power_profile')
        else
            power_profile = define_filter_settings(data_info,data_obj.ROI,data_obj.data,profile_range_mm,profile_step_mm);
            lowpass_cutoffs = lowpass_cutoffs(lowpass_cutoffs<(data_obj.filter_parameters.highpass*.9));
            filtersPwNumber = find_lowpassPwNumber(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
            save(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumber','power_profile')
        end
    end
    %% plot filter parameter
    
    f = figure();
    t = tiledlayout(1,3);
    f.Position = [100 100 1800 400];
    
    nexttile;
    plot(power_profile.scale_mm,power_profile.values,'DisplayName','power profile unfiltered map');
    hold on
    plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','lowpass cutoff')
    hold on
    plot([data_obj.filter_parameters.highpass data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','highpass cutoff')

    xlabel('Scale in mm')
    ylabel('Power')
    xlim([0.01 2])
    ylim([min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
    set(gca,'fontsize',15)
    legend('Location','northwest')
    
    nexttile;
    
    plot(average_spacing_mm./power_profile.scale_mm,power_profile.values,'DisplayName','power profile unfiltered map');
    hold on
    plot([average_spacing_mm./data_obj.filter_parameters.lowpass average_spacing_mm./data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','lowpass cutoff')
    hold on
    plot([average_spacing_mm./data_obj.filter_parameters.highpass average_spacing_mm./data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','highpass cutoff')

    xlabel('Wavevector (1/Λ)')
    ylabel('Power')
    xlim([average_spacing_mm/2 average_spacing_mm./(data_obj.filter_parameters.lowpass*0.8)])
    ylim([min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
    set(gca,'fontsize',15)
    legend()
    
    %Wavevector (1/Λ)

    
    nexttile;
    
    plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2),'DisplayName','pinwheel data')
    hold on
    plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all')],'DisplayName','lowpass cutoff')

    xlabel('Scale in mm')
    ylabel('# pinwheels')
    set(gca,'fontsize',15)
    legend()
    
    FigFile = [folder 'tFilterParameter' data_obj.info.ID];
    print(f,'-depsc', [FigFile '.eps'])
    savefig(f,[FigFile '.fig'])
    
end
