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
    [PwInfo.NumberPw,PwInfo.aniso,PwInfo.x_angle,PwInfo.PWxList,PwInfo.PWyList,PwInfo.signList, PwInfo.contours] = find_pinwheels(z,0,data_obj.ROI);
    PwInfo.average_spacing_mm = mean(mean(local_spacing_mm(oldROI == 1)));
    PwInfo.NumHypercolumns = sum(data_obj.ROI,'all')/(data_obj.info.pix_per_mm*PwInfo.average_spacing_mm)^2;
    PwInfo.MeanPwDensity = PwInfo.NumberPw/PwInfo.NumHypercolumns;
    
    %% get Local PwDensity Fixed Filter
    sigma = 0.1;
    PwInfo.LocalPwDensityFixedFilter=getLocalPwDensityFixedFilter(data_obj,PwInfo,local_spacing_mm,sigma);
    PwInfo.WeightedPwDensityFixedFilter = mean(PwInfo.LocalPwDensityFixedFilter(data_obj.ROI));
    
end

function LocalPwDensityFixedFilter = getLocalPwDensityFixedFilter(data_obj,PwInfo,local_spacing_mm,sigma)
    disp('calc LocalPwDensity with FixedFilter')
    average_spacing_mm = mean(local_spacing_mm(data_obj.ROI));
    local_pw_dens = put_gaussians(size(data_obj.ROI,1),size(data_obj.ROI,2), PwInfo.PWxList, PwInfo.PWyList,average_spacing_mm*data_obj.info.pix_per_mm,sigma,data_obj.ROI);
    local_pw_dens = local_pw_dens./sum(local_pw_dens(data_obj.ROI)).*PwInfo.NumberPw;
    LocalPwDensityFixedFilter = local_pw_dens.*(local_spacing_mm*data_obj.info.pix_per_mm).^2;
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