
function PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta)
    
    if nargin <6
        do_plotting=0;
    end
    if nargin <7
        llp_cutoffs = linspace(0.01, 1,100);
    end
    if nargin <8
        beta=0.5;
    end
    
    PwInfoFile = [DataFolder 'PwInfo_' data_obj.info.ID '.mat'];
    if isfile(PwInfoFile)
        load(PwInfoFile,'PwInfo')
    else
        
        data_obj.set_ROI(newROI)
        
        local_w = local_spacing_mm ;%* measure;
        local_w(newROI == 0) = 0;
        average_w = mean(mean(local_w(newROI == 1)));

        disp('Starting to estimate local pinwheel density with different cutoffs... ');
        PwInfo = estimate_local_pw_densityManuel(data_obj,average_w,local_w,llp_cutoffs,beta, do_plotting);

        disp('Analyzing pw NN distance statistics....');
        [PwInfo.d, PwInfo.d_eq, PwInfo.d_op] = compute_nn_pw_distances(PwInfo.PWxList,PwInfo.PWyList,PwInfo.signlist);


        PwInfo.map_area = sum(data_obj.ROI(:))/(average_w*data_obj.info.pix_per_mm)^2;


        % Now Compute count variance of pw numbers in subregions of different sizes
        % using the weighted pinwheel position matrix computed in
        % estimate_local_pw_density.m

        % Generate a weighted pw_pos_matrix for all pinwheels first
        PwInfo.pw_pos_matrix = PwInfo.pw_pos_matrix_plus + PwInfo.pw_pos_matrix_minus;
        weights = (average_w./local_w);
        PwInfo.weighted_pw_pos_matrix = PwInfo.pw_pos_matrix .*weights.*weights;

        % Compute count variance
        ITER = 500;
        [PwInfo.circ_areas, PwInfo.n] = gimme_nv_roi_with_local_pw_dens(PwInfo.weighted_pw_pos_matrix, ITER,average_w, data_obj.ROI);

        %% get Pinwheel pos stats
        tracker_obj = pinwheel_tracker;
        simple_track=true;
        [PwInfo.pinwheel_stats,PwInfo.pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
        
        %% Save PwInfo
        save(PwInfoFile,'PwInfo')
        
    end
    
    %% get PwDensityCIs
    alpha = 0.05;
    if getCI == true
        CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
        if isfile(CIPwFile)
            load(CIPwFile,'CI_PwDensities','CI_local_PwDensities','alpha')
            disp(['loaded data for alpha= ' num2str(alpha)])
            
        else
            local_w = local_spacing_mm ;
            local_w(newROI == 0) = 0;
            average_w = mean(mean(local_w(newROI == 1)));
            
            num_boot_samples = size(data_obj.samples_array,3);
            %% get PwDensits of bootstraped map 
            bootstat_PwDensities = zeros(1,size(data_obj.samples_array,3)-1);
            bootstat_local_PwDensities = zeros(sum(data_obj.ROI(:)),num_boot_samples-1);
            for ii = 2:num_boot_samples
                PwInfo.CI_Bootstrap{ii-1} = estimate_local_pw_densityManuel(data_obj,average_w,local_w,llp_cutoffs,beta, false,ii);
                bootstat_PwDensities(ii-1) = PwInfo.CI_Bootstrap{ii-1}.pw_dens;              
                bootstat_local_PwDensities(:,ii-1) = data_obj.array2vector(PwInfo.CI_Bootstrap{ii-1}.automated_pw_density);
                
            end
            
            %% get PwDensits of jackknife samples
            samples_array = data_obj.samples_array;
            data_obj.prepare_jackknife_samples;
            jackstat_PwDensities = zeros(1,size(data_obj.samples_array,3));
            jackstat_local_PwDensities = zeros(sum(data_obj.ROI(:)),size(data_obj.samples_array,3));
            for ii=1:data_obj.data_parameters.num_blocks
                PwInfo.CI_Jackstat{ii} = estimate_local_pw_densityManuel(data_obj,average_w,local_w,llp_cutoffs,beta, false,ii);
                jackstat_PwDensities(ii) = PwInfo.CI_Jackstat{ii}.pw_dens;              
                jackstat_local_PwDensities(:,ii) = data_obj.array2vector(PwInfo.CI_Jackstat{ii}.automated_pw_density);
            end
            data_obj.set_samples_array(samples_array);
            
            %% Calc CIs
            CI_PwDensities = bootstrap_ci(bootstat_PwDensities,PwInfo.pw_dens,jackstat_PwDensities,alpha);
            
            CI_local_PwDensitiesVector = bootstrap_ci(bootstat_local_PwDensities,data_obj.array2vector(PwInfo.automated_pw_density),jackstat_local_PwDensities,alpha);
            CI_local_PwDensities = zeros([size(data_obj.ROI) 2]);
            CI_local_PwDensities(:,:,1) = data_obj.vector2array(CI_local_PwDensitiesVector(:,1));
            CI_local_PwDensities(:,:,2) = data_obj.vector2array(CI_local_PwDensitiesVector(:,2));
            
            save(CIPwFile,'CI_PwDensities','CI_local_PwDensities','PwInfo','alpha')
        end
        PwInfo.CI_PwDensities=CI_PwDensities;
        PwInfo.CI_local_PwDensities = CI_local_PwDensities;
        PwInfo.alpha = alpha;
    end

end