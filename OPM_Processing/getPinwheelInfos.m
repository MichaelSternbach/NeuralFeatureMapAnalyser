function PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta)
    
    if nargin <6
        do_plotting=0;
    end
    if nargin <7
        average_spacing_mm = mean(mean(local_spacing_mm(newROI == 1)));
        llp_cutoffs = linspace(0.2*average_spacing_mm, 1*average_spacing_mm,50);
    end
    if nargin <8
        beta=0.5;
    end
    
    PwInfoFile = [DataFolder 'PwInfo_' data_obj.info.ID '.mat'];
    if isfile(PwInfoFile)
        disp('load Pw spacial info')
        load(PwInfoFile,'PwInfo')
    else
        
        %% get Pw spacial info
        disp('get Pw spacial info')
        PwInfo = getSpacialInfoPw(data_obj,local_spacing_mm,newROI,do_plotting,llp_cutoffs,beta);
        
        %% get Pinwheel pos stats
        tracker_obj = pinwheel_tracker;
        simple_track=true;
        [PwInfo.pinwheel_stats,PwInfo.pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
        
        
        
        
        %% Save PwInfo
        save(PwInfoFile,'PwInfo')
        
    end
    
    %% get PwDensityCIs
   disp('getCI')
   disp(getCI)
    %alpha = 0.05;
    if getCI == 1
        CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
        if isfile(CIPwFile)
            disp('load PwDensityCIs')
            load(CIPwFile,'CI_PwDensities','alpha')
            disp(['loaded data for alpha= ' num2str(alpha)])
            
        else
            disp('get PwDensityCIs')
            
            %% load spacing data
            CISpacingFile = [DataFolder 'CI_MapSpacing_' data_obj.info.ID '.mat'];
            load(CISpacingFile,'local_spacings_mm','local_spacingsJS_mm','newROIsBS','newROIsJS','alpha')
            disp(['loaded data for alpha= ' num2str(alpha)])
      
            
            num_boot_samples = size(data_obj.samples_array,3);
            %% get PwDensits of bootstraped map 
            disp('get PwDensits of bootstraped map ')
            PwInfosBS{1} = PwInfo;
            for ii = 2:num_boot_samples
                disp(['BS' num2str(ii)])
                PwInfosBS{ii} = getSpacialInfoPw(data_obj,local_spacings_mm{ii},newROIsBS{ii},false,llp_cutoffs,beta,ii);                
            end
            
            %% get PwDensits of jackknife samples
            disp('get PwDensits of jackknife samples')
            samples_array = data_obj.samples_array;
            data_obj.prepare_jackknife_samples;
            for ii=1:data_obj.data_parameters.num_blocks
                disp(['JS' num2str(ii)])
                PwInfosJS{ii} = getSpacialInfoPw(data_obj,local_spacingsJS_mm{ii},newROIsJS{ii},false,llp_cutoffs,beta,ii); 
            end
            data_obj.set_samples_array(samples_array);
            
            %% Calc CIs
            
            CI_PwDensities.PwDensityPosEstimate = getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'PwDensityPosEstimate',false,alpha);
            CI_PwDensities.PwDensityPlateuFit = getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'PwDensityPlateuFit',false,alpha);
            CI_PwDensities.LocalPwDensityPlateuFit = getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,'automated_pw_density',true,alpha);
                        
            save(CIPwFile,'CI_PwDensities','PwInfosBS','PwInfosJS','alpha')
        end
        PwInfo.CI_PwDensities=CI_PwDensities;
        PwInfo.alpha = alpha;
    end

end

