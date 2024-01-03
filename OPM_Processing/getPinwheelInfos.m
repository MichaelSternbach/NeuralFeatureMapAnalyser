
function PwInfo= getPinwheelInfos(data_obj,local_spacing_mm,DataFolder,newROI,getCI,do_plotting,llp_cutoffs,beta)
    
    if nargin <6
        do_plotting=0;
    end
    if nargin <7
        average_spacing_mm = mean(mean(local_spacing_mm(data_obj.ROI == 1)));
        llp_cutoffs = linspace(0.2*average_spacing_mm, 1*average_spacing_mm,50);
    end
    if nargin <8
        beta=0.5;
    end
    
    PwInfoFile = [DataFolder 'PwInfo_' data_obj.info.ID '.mat'];
    if isfile(PwInfoFile)
        load(PwInfoFile,'PwInfo')
    else
        
        %% get Pw spacial info
        
        PwInfo = getSpacialInfoPw(data_obj,local_spacing_mm,newROI,do_plotting,llp_cutoffs,beta);
        
        %% get Pinwheel pos stats
        tracker_obj = pinwheel_tracker;
        simple_track=true;
        [PwInfo.pinwheel_stats,PwInfo.pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);

        %% Save PwInfo
        save(PwInfoFile,'PwInfo')
        
    end
    
    %% get PwDensityCIs
    %alpha = 0.05;
    if getCI == true
        CIPwFile = [DataFolder 'CI_PwDensity_' data_obj.info.ID '.mat'];
        if isfile(CIPwFile)
            load(CIPwFile,'CI_PwDensities','alpha')
            disp(['loaded data for alpha= ' num2str(alpha)])
            
        else
            
            %% load spacing data
            CISpacingFile = [DataFolder 'CI_MapSpacing_' data_obj.info.ID '.mat'];
            load(CISpacingFile,'local_spacings_mm','local_spacingsJS_mm','newROIsBS','newROIsJS','alpha')
            disp(['loaded data for alpha= ' num2str(alpha)])
      
            
            num_boot_samples = size(data_obj.samples_array,3);
            %% get PwDensits of bootstraped map 
            PwInfosBS{1} = PwInfo;
            for ii = 2:num_boot_samples
                PwInfosBS{ii} = getSpacialInfoPw(data_obj,local_spacings_mm{ii},newROIsBS{ii},false,llp_cutoffs,beta,ii);                
            end
            
            %% get PwDensits of jackknife samples
            samples_array = data_obj.samples_array;
            data_obj.prepare_jackknife_samples;
            for ii=1:data_obj.data_parameters.num_blocks
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



function PwInfo = getSpacialInfoPw(data_obj,local_spacing_mm,newROI,do_plotting,llp_cutoffs,beta,sample)

    if nargin < 7
        sample=1;
    end

    oldROI = data_obj.ROI;

    %% calc local Pw Density and other Pw position stats

    PwInfo = calcLocalPwDensityAndPosStats(data_obj,local_spacing_mm,newROI,do_plotting,llp_cutoffs,beta,sample);


    %% get mean pinwheel density full ROI
    data_obj.set_ROI(oldROI)
    z = data_obj.filter_map(data_obj.read_map(sample));
    [PwInfo.NumberPw,~,~,~,~,~, ~] = find_pinwheels(z,0,data_obj.ROI);
    PwInfo.average_spacing_mm = mean(mean(local_spacing_mm(oldROI == 1)));
    PwInfo.NumHypercolumns = sum(data_obj.ROI,'all')/(data_obj.info.pix_per_mm*PwInfo.average_spacing_mm)^2;
    PwInfo.MeanPwDensity = PwInfo.NumberPw/PwInfo.NumHypercolumns;

end

function PwInfo = calcLocalPwDensityAndPosStats(data_obj,local_spacing_mm,newROI,do_plotting,llp_cutoffs,beta,sample)

    if nargin < 7
        sample=1;
    end

    %% calculate local pw density
    data_obj.set_ROI(newROI)

    local_w = local_spacing_mm ;%* measure;
    local_w(newROI == 0) = 0;
    average_w = mean(mean(local_w(newROI == 1)));

    disp('Starting to estimate local pinwheel density with different cutoffs... ');
    PwInfo = estimate_local_pw_densityManuel(data_obj,average_w,local_w,llp_cutoffs,beta, do_plotting, sample);

    disp('Analyzing pw NN distance statistics....');
    try
        [PwInfo.d, PwInfo.d_eq, PwInfo.d_op] = compute_nn_pw_distances(PwInfo.PWxList,PwInfo.PWyList,PwInfo.signlist);
    catch
        disp('ERROR in compute_nn_pw_distances!')
    end

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

end