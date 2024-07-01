function AnimalAllFiltering(animal,experiment_Num,AnimalDataFolder,DataFolderMain,FigureFolder)
    close all
    
    %% parameter
    experiment_Num = checkFormatNum(experiment_Num);
    
    
    rectangle = false;
    
    
    %% figure folder
    mkdir(FigureFolder)
    

    FigureFile1 = [FigureFolder animal '_PwPlateau'];
    FigureFile2 = [FigureFolder animal '_PowerSpectrum'];
    FigureFile3 = [FigureFolder animal '_PwDensPlateau'];
    
    rm_cmd = ['rm -f ' FigureFile1 '.ps'];
    system(rm_cmd)
    rm_cmd = ['rm -f ' FigureFile2 '.ps'];
    system(rm_cmd)
    
    
    
    %% make Figure
% 
%     f1 = figure(1);
%     t1 = tiledlayout(3,3);
%     f1.Position = [100 100 1000 8000];
    
    f2 = figure(2);
    t2 = tiledlayout(3,3);
    f2.Position = [100 100 1000 8000];
    
    f3 = figure(3);
    t3 = tiledlayout(3,3);
    f3.Position = [100 100 1000 8000];

    for experiment_num = 1: experiment_Num
        
        %% data folder
        DataFolder = [DataFolderMain lower(animal) '/' lower(animal) num2str(experiment_num) '/'];
        
        %% animal
        [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,1,false,AnimalDataFolder);
        [average_spacing_mm,local_spacing_mm,~] =  getColumnsSpacing(data_obj,DataFolder,false);



        %% load pinwheel data

        llp_cutoffs = linspace(0.01, 1,100);
        beta=0.5;

        z = data_obj.filter_map(data_obj.read_map());
        %[~,~,~,PWxList,PWyList,~, ~] = find_pinwheels(z,0,data_obj.ROI);





        %% pw Plateau
        DataAndFilterFile = [DataFolder 'FilterFile.mat'];

%         if isfile(DataAndFilterFile) %&& false 
%             load(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumber')
%         else
% 
%             lowpass_cutoffs = linspace(0.2*average_spacing_mm,average_spacing_mm*0.7,100);
%             filtersPwNumber = find_lowpassPwNumber(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
%             save(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumber')
%         end
        
        %% Filter data
        DataAndFilterFile = [DataFolder 'FilterFile.mat'];

        if isfile(DataAndFilterFile) && false 
            load(DataAndFilterFile,'lowpass_cutoffs','filtersPwNumber','power_profile')
        else
            lowpass_cutoffs = linspace(0.1,0.6,50);%linspace(0.2*average_spacing_mm,average_spacing_mm*.8,200);
            profile_range_mm = [0.01 2]; 
            profile_step_mm = 0.01;
            %disp(lowpass_cutoffs)
            power_profile = define_filter_settings(data_info,data_obj.ROI,data_obj.data,profile_range_mm,profile_step_mm);
            lowpass_cutoffs = lowpass_cutoffs(lowpass_cutoffs<(data_obj.filter_parameters.highpass*.9));
%             filtersPwNumber = find_lowpassPwNumber(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
            [filtersPwDens,pltPwDens, paramsPwDens, sselPwDens] = find_lowpassPwNumberFit(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs,true,local_spacing_mm);
            %[filtersPwNumber,plt, params, ssel] = find_lowpassPwNumberFit(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
            %save(DataAndFilterFile,'data_obj','data_info','lowpass_cutoffs','filtersPwNumber','power_profile')
        end
        
%         yy = PlateauFit(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),plt, params, ssel);
%         figure(1);
%         nexttile;
%         plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2),'x','DisplayName','Pinwheel Data')
%         hold on;
%         plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),yy,'r','DisplayName','fit')
%         hold on
%         minPw = min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2))*0.9;
%         maxPw = max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2))*1.1;
%         plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[minPw maxPw],'b','DisplayName','Lowpass Cutoff')
%         %min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all')
%         %xlim([0.2 .8].*average_spacing_mm)
%         ylim([minPw maxPw])
%         xlim([0.1 0.4])
%         title(data_obj.info.ID)
%         legend()
%         xlabel('Scale [mm]')
%         ylabel('# Pinwheels')
%         axis square;    
        
        
        yyPwDens = PlateauFit(filtersPwDens.global_plateau.lowpass_vs_density(:,1),pltPwDens, paramsPwDens, sselPwDens);
        figure(3);
        nexttile;
        plot(filtersPwDens.global_plateau.lowpass_vs_density(:,1),filtersPwDens.global_plateau.lowpass_vs_density(:,2),'x','DisplayName','Pinwheel Data')
        hold on;
        plot(filtersPwDens.global_plateau.lowpass_vs_density(:,1),yyPwDens,'r','DisplayName','fit')
        hold on
        minPw = min(filtersPwDens.global_plateau.lowpass_vs_density(:,2))*0.9;
        maxPw = max(filtersPwDens.global_plateau.lowpass_vs_density(:,2))*1.1;
        plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[minPw maxPw],'b','DisplayName','Lowpass Cutoff')
        %min(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all') max(filtersPwNumber.global_plateau.lowpass_vs_density(:,2),[],'all')
        %xlim([0.2 .8].*average_spacing_mm)
        ylim([minPw maxPw])
        xlim([0.1 0.4])
        title(data_obj.info.ID)
        legend()
        xlabel('Scale [mm]')
        ylabel('Pinwheel Density')
        axis square;    
        
        %% powerspectrum OPM
        figure(2);
        nexttile;

%         plot(average_spacing_mm./power_profile.scale_mm,power_profile.values,'DisplayName','Power Profile Unfiltered Map');
%         hold on
%         plot([average_spacing_mm./data_obj.filter_parameters.lowpass average_spacing_mm./data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','Lowpass Cutoff')
%         hold on
%         plot([average_spacing_mm./data_obj.filter_parameters.highpass average_spacing_mm./data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','Highpass Cutoff')

        plot(1./power_profile.scale_mm,power_profile.values,'DisplayName','Power Profile Unfiltered Map');
        hold on
        plot([1./data_obj.filter_parameters.lowpass 1./data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','Lowpass Cutoff')
        hold on
        plot([1./data_obj.filter_parameters.highpass 1./data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','Highpass Cutoff')


        xlabel('Wavevector [1/mm]')
        ylabel('Power')
        xlim([1 1./(data_obj.filter_parameters.lowpass*0.8)])
        ylim([0 250])% max(power_profile.values,[],'all')])
        %set(gca,'fontsize',15)
        axis square
        title(data_obj.info.ID)
        %axis tight
        legend('Location','northoutside')
      
      
        %% save figure
        %print(f,'-depsc', [FigureFile '.eps'])
        %savefig(f,[FigureFile '.fig'])
        
    end
    

   % print(f1, '-dpsc','-fillpage', [FigureFile1 '.ps'])%'-append',
    print(f2, '-dpsc','-fillpage', [FigureFile2 '.ps'])%'-append',
    print(f3, '-dpsc','-fillpage', [FigureFile3 '.ps'])%'-append',
    
    %% finished
    disp('saved & finished')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return 
end

function yy = PlateauFit(x,plt, params, ssel)
    %plt = [mnplt_min mxplt_min]
    mnplt_min = plt(1);
    mxplt_min = plt(2);
    % Compute the actual fit function
    if isequal(ssel, [1, 1, 1])
        yy = params(1) + (-x+mnplt_min).*((-x+mnplt_min)>0)*params(2) + (-x+mxplt_min).*((-x+mxplt_min)<0)*params(3);       
    elseif isequal(ssel, [1, 0, 1])
        yy = params(1) +  (-x+mxplt_min).*((-x+mxplt_min)<0)*params(2);
    elseif isequal(ssel, [1, 1, 0])
        yy = params(1) + (-x+mnplt_min).*((-x+mnplt_min)>0)*params(2);
    end

end



