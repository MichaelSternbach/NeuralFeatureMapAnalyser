function data_info = getFilterSettings(data_obj,data_info,folder,splitROI,resetFilter, lowpass_cutoffs,profile_range_mm)
    
    if nargin <4
        splitROI = false;
    end
    if ~isempty(splitROI) && (~islogical(splitROI) || (splitROI ~= false))
        if isnumeric(splitROI)
            NumSplitROI = splitROI;
            [splitROI,partialROIs] = getSplitROI(data_obj.ROI,NumSplitROI);
        else
            partialROIs = getPartialROIs(data_obj.ROI,splitROI);
        end
    elseif isempty(splitROI) || (islogical(splitROI) && ~splitROI)
        splitROI=cell(1);
        splitROI{1}.x = [1 size(data_obj.ROI,1)];
        splitROI{1}.y = [1 size(data_obj.ROI,2)];
        ROI=data_obj.ROI;
        partialROIs = reshape(ROI,[size(ROI) 1]);
    end
    
    if nargin <5
        resetFilter = false;
    end
    if nargin <6
        lowpass_cutoffs = 0.1:0.1:2;
    end
    
    if nargin <7
        profile_range_mm = 0.1:0.01: 2;
    end
%     if nargin <8 
%         profile_step_mm = 0.01;
%     end
    
    %% define file names
    DataAndFilterFile = [folder data_info.ID '.mat'];
    
    %% save ROI
    ROI = data_obj.ROI;
    %% calc x/y range
    [x_range,y_range] = getRangeXY_ROI(ROI);

    if length(splitROI) < 7
        color_list = {'r','g','b','c','m','y','k'};
    else
        %color_list = jet(length(splitROI));
        color_list = getColorList(length(splitROI));
    end
    
    %% input filter parameter
    if resetFilter || ~isfield(data_info,'settings') || ~isfield(data_info.settings,'lowpass_mm') || ~isfield(data_info.settings,'highpass_mm')
        
        %% calculate power spectrum
        power_profile = define_filter_settings(data_info,data_obj.ROI,data_obj.data,profile_range_mm);

        figure
        plot(power_profile.scale_mm,power_profile.values);
        xlabel('Scale in mm')
        ylabel('Power')
        xlim([0 max(profile_range_mm)])
        set(gca,'fontsize',15)
        
        
        %% input filter settings
        prompt = "What is the lowpass filter value in mm? ";
        lowpass_mm = input(prompt);
        data_obj.set_filter_parameters('lowpass',lowpass_mm)
        data_info.settings.lowpass_mm = lowpass_mm;

        prompt = "What is the highpass filter value in mm? ";
        highpass_mm = input(prompt);
        data_obj.set_filter_parameters('highpass',highpass_mm)
        data_info.settings.highpass_mm = highpass_mm;

        close all
        
%         %% plot OPM 
% 
%         f0 = figure;
%         z1 = data_obj.filter_map(data_obj.read_map()); %map
%         plot_map(z1,data_obj.ROI,0,1)

        %% determine lowpass filter from pinwheel plateau

        lowpass_cutoffs = lowpass_cutoffs(lowpass_cutoffs<(highpass_mm*.9));
        f1=figure;

        filtersPwNumbers = cell(length(splitROI),1);
        for ii = 1:length(splitROI)
            data_obj.set_ROI(partialROIs(:,:,ii));
            filtersPwNumbers{ii} = find_lowpassPwNumber(data_obj,data_info,highpass_mm,lowpass_mm,lowpass_cutoffs);
            plot(filtersPwNumbers{ii}.global_plateau.lowpass_vs_density(:,1),filtersPwNumbers{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)],'Color',color_list{ii})
            hold on
        end
        % filtersPwNumber = find_lowpassPwNumber(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
        % plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2),'DisplayName','pinwheel data')
        
        legend()
        xlabel('lowpassfilter cutoff [mm]')
        ylabel('pinwheel number')

        
        f2 = figure;
        %% plot map
        plot_map(data_obj.filter_map(data_obj.read_map()),ROI,0,1)
        %% plot partial ROIs (rectangles)
        for ii = 1:length(splitROI)
            rectangle('Position',[splitROI{ii}.x(1),splitROI{ii}.y(1),splitROI{ii}.x(2)-splitROI{ii}.x(1),splitROI{ii}.y(2)-splitROI{ii}.y(1)],'EdgeColor',color_list{ii},'LineWidth',2)
        end
        xlim([x_range(1) x_range(2)])
        ylim([y_range(1) y_range(2)])

        prompt = "What is the lowpass filter value in mm based on the pinwheel density plateau? ";
        lowpass_mm = input(prompt);
        data_obj.set_filter_parameters('lowpass',lowpass_mm)
        data_info.settings.lowpass_mm = lowpass_mm;
        close(f1)
        close(f2)

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
        
        
        % figure 
        % plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2),'DisplayName','pinwheel data')
        % hold on
        % plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all')],'DisplayName','lowpass cutoff')

        % xlabel('Scale in mm')
        % ylabel('# pinwheels')
        % set(gca,'fontsize',15)
        % legend()

        figure;
        for ii = 1:length(splitROI)
            plot(filtersPwNumbers{ii}.global_plateau.lowpass_vs_density(:,1),filtersPwNumbers{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)],'Color',color_list{ii})
            hold on
        end
        % plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all')],'DisplayName','lowpass cutoff')

        legend()
        xlabel('lowpassfilter cutoff [mm]')
        ylabel('pinwheel number')


        figure;
        %% plot map
        plot_map(data_obj.filter_map(data_obj.read_map()),ROI,0,1)
        %% plot partial ROIs (rectangles)
        for ii = 1:length(splitROI)
            rectangle('Position',[splitROI{ii}.x(1),splitROI{ii}.y(1),splitROI{ii}.x(2)-splitROI{ii}.x(1),splitROI{ii}.y(2)-splitROI{ii}.y(1)],'EdgeColor',color_list{ii},'LineWidth',2)
        end
        xlim([x_range(1) x_range(2)])
        ylim([y_range(1) y_range(2)])


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
        save(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumbers','power_profile','partialROIs','splitROI')
    else
        %% use existing filter parameter and data
        if isfile(DataAndFilterFile)
            load(DataAndFilterFile,'filtersPwNumbers','power_profile','partialROIs','splitROI')
        else
            power_profile = define_filter_settings(data_info,data_obj.ROI,data_obj.data,profile_range_mm);
            lowpass_cutoffs = lowpass_cutoffs(lowpass_cutoffs<(data_obj.filter_parameters.highpass*.9));
            %filtersPwNumber = find_lowpassPwNumber(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
            filtersPwNumbers = cell(length(splitROI),1);
            for ii = 1:length(splitROI)
                data_obj.set_ROI(partialROIs(:,:,ii));
                filtersPwNumbers{ii} = find_lowpassPwNumber(data_obj,data_info,data_info.settings.highpass_mm,data_info.settings.lowpass_mm,lowpass_cutoffs);
            end
            %save(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumbers','power_profile','partialROIs','splitROI')
        end
    end
    %% plot filter parameter
    
    f = figure();
    t = tiledlayout(1,4);
    f.Position = [100 100 1800 400];
    title(t,'Filter Choice')
    
    %% plot power profile
    nexttile;
    plot(power_profile.scale_mm,power_profile.values,'DisplayName','power profile unfiltered map');
    hold on
    plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','lowpass cutoff')
    hold on
    plot([data_obj.filter_parameters.highpass data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','highpass cutoff')

    xlabel('Scale in mm')
    ylabel('Power')
%     xlim([0.01 2])
%     ylim([min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
    set(gca,'fontsize',15)
    %legend('Location','northwest')
    

    %% plot power spectrum
    nexttile;
    
    plot(1./power_profile.scale_mm,power_profile.values_kspace,'DisplayName','power spectrum unfiltered map');
    hold on
    plot([1./data_obj.filter_parameters.lowpass 1./data_obj.filter_parameters.lowpass],[min(power_profile.values_kspace,[],'all') max(power_profile.values_kspace,[],'all')],'DisplayName','lowpass cutoff')
    hold on
    plot([1./data_obj.filter_parameters.highpass 1./data_obj.filter_parameters.highpass],[min(power_profile.values_kspace,[],'all') max(power_profile.values_kspace,[],'all')],'DisplayName','highpass cutoff')

    xlabel('Wavevector (1/mm)')
    ylabel('Power')
%     xlim([1./(data_obj.filter_parameters.highpass*1.3) 1./(data_obj.filter_parameters.lowpass*0.8)])
%     ylim([min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
    set(gca,'fontsize',15)
    %legend()

    %% plot pw plateau
    nexttile;

    for ii = 1:size(partialROIs,3)
        plot(filtersPwNumbers{ii}.global_plateau.lowpass_vs_density(:,1),filtersPwNumbers{ii}.global_plateau.lowpass_vs_density(:,2),'DisplayName',['pinwheel data ROI ' num2str(ii)],'Color',color_list{ii})
        hold on
    end
    plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(filtersPwNumbers{1}.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumbers{1}.global_plateau.lowpass_vs_density(:,2),[],'all')],'DisplayName','lowpass cutoff')

    %legend()
    xlabel('lowpassfilter cutoff [mm]')
    ylabel('pinwheel number')


    nexttile;
    %% plot map
    data_obj.set_ROI(ROI);
    plot_map(data_obj.filter_map(data_obj.read_map()),ROI,0,1)
    hold on
    %% plot partial ROIs (rectangles)
    for ii = 1:length(splitROI)
        rectangle('Position',[splitROI{ii}.x(1),splitROI{ii}.y(1),splitROI{ii}.x(2)-splitROI{ii}.x(1),splitROI{ii}.y(2)-splitROI{ii}.y(1)],'EdgeColor',color_list{ii},'LineWidth',2)
    end
    xlim([x_range(1) x_range(2)])
    ylim([y_range(1) y_range(2)])

    %% save figure
    
    FigFile = [folder 'FilterParameter' data_obj.info.ID];
    print(f,'-depsc', [FigFile '.eps'])
    savefig(f,[FigFile '.fig'])
    
end

%% functions 

function [splitROI,partialROIs] = getSplitROI(ROI,NumSplitROI)
    %% get ranges ROI
    [x_range,y_range] = getRangeXY_ROI(ROI);
    %% get split ROI
    % splitROI = cell(NumSplitROI*NumSplitROI,1);
    % partialROIs = zeros(size(ROI,1),size(ROI,2),NumSplitROI*NumSplitROI);
    x_split = round(linspace(x_range(1),x_range(2),NumSplitROI+1));
    y_split = round(linspace(y_range(1),y_range(2),NumSplitROI+1));
    size_splitROI = (x_split(2)-x_split(1))*(y_split(2)-y_split(1));
    n=1;
    for ii = 1:NumSplitROI
        for jj = 1:NumSplitROI
            % %% get splitROI
            % splitROI{ii}.x = [x_split(ii) x_split(ii+1)];
            % splitROI{ii}.y = [y_split(jj) y_split(jj+1)];
    
            % %% calc ROI
            % partialROIs(:,:,ii) = getPartialROI(splitROI{ii}.x,splitROI{ii}.y,ROI);
            % %% get splitROI
            % splitROI{(ii-1)*NumSplitROI+jj}.x = [x_split(ii) x_split(ii+1)];
            % splitROI{(ii-1)*NumSplitROI+jj}.y = [y_split(jj) y_split(jj+1)];

            % %% calc ROI
            % partialROIs(:,:,((ii-1)*NumSplitROI+jj)) = getPartialROI(splitROI{(ii-1)*NumSplitROI+jj}.x,splitROI{(ii-1)*NumSplitROI+jj}.y,ROI);


            %% get splitROI
            splitROI_x = [x_split(ii) x_split(ii+1)];
            splitROI_y = [y_split(jj) y_split(jj+1)];

            %% calc ROI
            partialROI = getPartialROI(splitROI_x,splitROI_y,ROI);


            %% check size of splitROI
            if (sum(partialROI,'all')/size_splitROI) > 0.8            
                splitROI{n}.x = splitROI_x;
                splitROI{n}.y = splitROI_y;
                partialROIs(:,:,n) = partialROI;
                n = n+1;
            end
        end
    end

end

function partROI = getPartialROI(x,y,ROI)
    %% get Meshrid size ROI
    [X,Y] = meshgrid(1:size(ROI,2),1:size(ROI,1));
    %partROI = (X >= x_split(ii) & X <= x_split(ii+1)) & (Y >= y_split(ii) & Y <= y_split(ii+1));
    partROI_bool = ((X >= x(1) & X <= x(2)) & (Y >= y(1) & Y <= y(2)));
    partROI = zeros(size(ROI));
    partROI(partROI_bool)=1;
    partROI = (partROI.*ROI);
end                                                                                                                                                                                                                

function partialROIs = getPartialROIs(ROI,splitROI)
    partialROIs = zeros(size(ROI,1),size(ROI,2),length(splitROI));
    for ii = 1:length(splitROI)
        partialROIs(:,:,ii) = getPartialROI(splitROI{ii}.x,splitROI{ii}.y,ROI);
    end
end

%% get colorlist with large difference bewtween colors as cell
function color_list = getColorList(n)
    color_list = cell(n,1);
    for ii = 1:n
        color_list{ii} = [rand rand rand];
        for jj = 1:ii-1
            while sum(abs(color_list{ii}-color_list{jj})) < 0.5
                color_list{ii} = [rand rand rand];
            end
        end
    end
end 