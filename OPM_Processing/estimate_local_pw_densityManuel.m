function map = estimate_local_pw_densityManuel(map,hhp, llp_cutoffs,beta,do_plotting)    
%
% FUNCTION CALLED BY ANALYZE_SINGLE_MAP.M
%
% USAGE : 
% map = estimate_local_pw_density(map,hhp, llp_cutoffs,beta,do_plotting)    
%
% This function implements the automated pinwheel density estimation method introduced in Kaschube et al., Science 2010.
% The map is first highpass filtered and then low-pass filtered with
% different cut-offs. Pinhweels are identified in each of the filtered maps
% pinwheel density is estimated by a piece-wise linear fit and identifying
% a plateau  (see Kaschube et al., Science 2010 Supporting information)
%
% Actual pinwheel positions are estimated using filter values at the center of the
% plateaus for each point and determine the existence of pinwheels arond those pixels for
% these filter values
%
%
%
% intermediate results are stored during function execution in the analysis
% dir
%
%
%
% INPUT PARAMETERS:
% map          ... data structure containing all information gather so far (local wavelength etc.) (see help ANALYZE_SINGLE_MAP.M for more info
%
%  llp         ... lowpass cutoff frequency in pixels
%  hhp         ... highpass cutoff frequency in pixels
%                   used subsequently for processing
%  beta         ... "Temperature" for the fermi-filter if used,
%                   beta is measured in units of k_hp or k_lp (see Kaschube et al., Science 2010)
%
%
%  do_plotting   ... flag, either zero or one, chose 1 for visual output of
%                   filter results
% OUTPUT PARAMETERS:
%       map ... data structure containing all pw density analysis information (see help ANALYZE_SINGLE_MAP.M for more info
%
%       THE FOLLOWING FIELDS ARE ADDED to the map structure
%                   pw_x_lp ... local pinwheel density as a function of lowpass cutoff 
%           pw_pos_estimate ... array used for pinwheel position estimate
%      pw_pos_estimate_plus ... array used for pinwheel position estimate
%                               only positive charge pinwheels     
%     pw_pos_estimate_minus ... array used for pinwheel position estimate
%                               only negative charge pinwheels     
%      automated_pw_density ... local pinwheel density after plateau
%                               fitting
%                   pw_dens ... average pinwheel density of map
%             pw_pos_matrix ... matrix used to estimate the pinwheel positions [480x655 double]
%        pw_pos_matrix_plus ... matrix used to estimate the positions for pinwheel of positive charge 
%       pw_pos_matrix_minus ... matrix used to estimate the positions for pinwheel of negative charge 
%                   PWxList ... linear array for x-coordinates of pinwheels
%                   PWyList ... linear array for y-coordinates of pinwheels
%                  signlist ... linear array for charge of pinwheels (either -1 or 1)
%
%
%  Copyright (C) 2014, Max-Planck-Institut for Dynamics and Self-organization, The  
%  Nonlinear Dynamics Group. This software may be used, copied, or 
%  redistributed as long as it is not sold and this copyright notice is 
%  reproduced on each copy made. This routine is provided as is without 
%  any express or implied warranties whatsoever.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if beta <= 0
        error('Choose a beta value that is larger than one. beta = 0.05 is recommended.');
        
    end
        
    data_tmp = map.z_shrunk; %% z_shrunk is raw map (only shrunk, without any filtering)

    %%%% Set everything outside of the ROI to zero
    data_tmp(map.new_roi == 0) = 0;
    
    
    %%%% Define arrays for pw position estimation    
    pw_x_lp = zeros(size(data_tmp,1),size(data_tmp,2),length(llp_cutoffs));

    pw_pos_estimate_plus = zeros(size(data_tmp,1),size(data_tmp,2),length(llp_cutoffs));
    pw_pos_estimate_minus = zeros(size(data_tmp,1),size(data_tmp,2),length(llp_cutoffs));


    %%%%%%%%%%% High pass filtering first, no low pass filtering
    [data_hp]= get_filtered_for_complex_fields(data_tmp,map.new_roi,'highpass',hhp,beta);
    data_hp(map.new_roi == 0) = 0;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Now filter  with different lowpass cutoffs
    if ~isdeployed && do_plotting
        scrsz = get(0,'ScreenSize');
        load colormap_min;
        figure(5)
        set(gcf, 'Position',[1 1 scrsz(3)/2 scrsz(4)/2]);
        figure(6)
        set(gcf, 'Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)/2]);
    end
    
    PWX_all_cutoffs = zeros(1000,length(llp_cutoffs));
    PWY_all_cutoffs = zeros(1000,length(llp_cutoffs));
    
    for ii = 1:length(llp_cutoffs)      

        data_tmp = data_hp;   

        % Lowpass Fermi filtering with different cut-offs
        [data_tmp]= get_filtered_for_complex_fields(data_tmp,map.new_roi, 'lowpass', llp_cutoffs(ii),beta);
        data_tmp(map.new_roi == 0) = 0;

        % now set variance to one
        data_tmp(map.new_roi == 1) = data_tmp(map.new_roi == 1)/std(data_tmp(map.new_roi == 1));


        %%% Finding pinwheels and computing local pw density
        [~,~,~,PWxList,PWyList,signList] = pw_finder_withsign_DEWRevision(data_tmp,1,map.new_roi,0);             
        
        PWX_all_cutoffs(1:length(PWxList),ii) = PWxList;
        PWY_all_cutoffs(1:length(PWyList),ii) = PWyList;
        
        %%% Creates array for pinwheel density estimation
        local_pw_dens = put_gaussians(size(data_tmp,1),size(data_tmp,2), PWxList, PWyList,map.average_w,1.0,map.new_roi);
        pw_x_lp(:,:,ii) = local_pw_dens.*map.local_w.^2;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This is for finding pinwheels
        % First make the pinwheel positions integer values to have them be
        % indices
        
        PWxList = floor(PWxList);
        PWyList = floor(PWyList);
        
        % Make sure that no pinwheel position is at the very border of
        % the image (these come out with sign = 0 from the pinwheel finder routine)
        PWxList = PWxList(PWyList > 0);
        signList = signList(PWyList > 0);
        PWyList = PWyList(PWyList > 0);
        
        PWyList = PWyList(PWxList > 0);
        signList = signList(PWxList > 0);
        PWxList = PWxList(PWxList > 0);
        
        
        PWxplus = PWxList(signList > 0);
        PWyplus = PWyList(signList > 0);
        iindex = sub2ind(size(data_tmp),PWyplus,PWxplus);
        tmp = zeros(size(data_tmp));
        tmp(iindex) = 1;
        pw_pos_estimate_plus(:,:, ii) = tmp;
        
        PWxminus = PWxList(signList < 0);
        PWyminus = PWyList(signList < 0);
        iindex = sub2ind(size(data_tmp),PWyminus,PWxminus);
        tmp = zeros(size(data_tmp));
        tmp(iindex) = 1;
        pw_pos_estimate_minus(:,:, ii) = tmp;

        if ~isdeployed
            figure(5)
            imagesc(angle(data_tmp));
            colormap(colormap_min);
            axis image;
            hold on;
            plot(PWxList,PWyList,'o','color','w','MarkerSize',5, 'Linewidth',2);
            title(['OPM and pw locations for for cut-off ' num2str(llp_cutoffs(ii)/map.measure) 'mm']);
            set(gca, 'ydir', 'normal', 'fontsize', 16);
            hold off;
            
            figure(6)
            imagesc(pw_x_lp(:,:,ii));
            colormap(jet);
            axis image;
            hold on;
            plot(PWxList,PWyList,'o','color','w','MarkerSize',5, 'Linewidth',2);
            title(['Local pw density for cut-off ' num2str(llp_cutoffs(ii)/map.measure) 'mm']);
            colorbar;   
            set(gca, 'ydir', 'normal', 'fontsize', 16);
            hold off;
            
        end

    end% For loop
    
    %%%% Just for storage! These are rather big matrices and do not need to
    %%%% be stored in map-structure
    intermediate_results.pw_x_lp = pw_x_lp;   
    
    intermediate_results.pw_pos_estimate_plus = pw_pos_estimate_plus;   
    intermediate_results.pw_pos_estimate_minus = pw_pos_estimate_minus;   
    
    intermediate_results.PWX_all_cutoffs = PWX_all_cutoffs;   
    intermediate_results.PWY_all_cutoffs = PWY_all_cutoffs;   

    
    clear x;
    clear y;
    clear xi;
    clear yi;
    clear pw_pos_estimate_plus;
    clear pw_pos_estimate_minus;

%     save([map.data_dirname 'analysis/map_analysis_data_intermediate_results.mat'], 'intermediate_results');
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%     
%     load([map.data_dirname 'analysis/map_analysis_data_intermediate_results.mat'], 'intermediate_results');
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Starting piecewise linear fitting of pinwheel plateaus...');   
 
    [XROI, YROI] = find(map.new_roi > 0);      
    automated_density = zeros(size(map.new_roi));
    
    pw_pos_matrix_plus = zeros(size(map.new_roi));
    pw_pos_matrix_minus = zeros(size(map.new_roi));
    
    progress = 0.09999;

    disp('0% done... ');
    
    for ii = 1:length(XROI)
        
        if ii/length(XROI) > progress
            disp([ num2str(round(progress*100)) '% done...']);
            progress = progress + .1;
        end
        y = zeros(1,length(llp_cutoffs));
        y(:) = intermediate_results.pw_x_lp(XROI(ii),YROI(ii),:);
        x = llp_cutoffs./map.local_w(XROI(ii),YROI(ii));
        
        %%% Only fit in the cut-off interval (0.2, 1)Lambda
        y = y((x > 0.2) & (x < 1));
        x = x((x > 0.2) & (x < 1));
        
        %%%% Piecewise linear fitting
        [pwd, plt] = fit_piecewise_linear_WK_rerevised(x,y,0);
        automated_density(XROI(ii),YROI(ii)) = pwd;    
        
        
        %%%% Now find center of plateau
        cent_plt = sum(plt)/2;
        
        % Find index of llp cut off closest to center of plateau
        cutoff = find(abs(x-cent_plt) == min(abs(x-cent_plt)),1);
        
        %%% assign value at plateau center to matrices used to find pinwheels
        %%% later
        pw_pos_matrix_plus(XROI(ii),YROI(ii)) = intermediate_results.pw_pos_estimate_plus(XROI(ii),YROI(ii),cutoff);
        pw_pos_matrix_minus(XROI(ii),YROI(ii)) = intermediate_results.pw_pos_estimate_minus(XROI(ii),YROI(ii),cutoff);
                                
        
    end
    
    
    % Weight the local pinwheel density, such that each hypercolumn gets
    % the same weight when averaging    
    weights = (map.average_w./map.local_w);
    weighted_pw_density = (automated_density.*weights.*weights);
    weighted_pw_density(automated_density<=0) = 0;
    lpwd = mean(weighted_pw_density(weighted_pw_density>0));

    
    % Store results in map structure
    map.automated_pw_density = automated_density;
    map.automated_weighted_pw_density = weighted_pw_density;
    map.pw_dens = lpwd;
    
    map.pw_pos_matrix_plus = pw_pos_matrix_plus;
    map.pw_pos_matrix_minus = pw_pos_matrix_minus;
    

%     % Intermediate save 
%     save([map.data_dirname 'analysis/map_analysis_data_new.mat'], 'map');
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     load([map.data_dirname 'analysis/map_analysis_data_new.mat'], 'map');

    
    % Clean up pinwheel matrices by eliminating pinwheels that are too
    % close to each other
   
    disp('Cleaning up pw positions ...');

    
    %%% Set minimal distance for threshold
    thd = 0.03*map.measure;   
    
    %%% Clean up positively charged pinwheels
    [indx indy] = find(map.pw_pos_matrix_plus > 0);
    ii = 1;
    while ii < length(indx)
        
        dists = sqrt((indx(ii) - indx).^2 + (indy(ii) - indy).^2);        
        ind = find(dists < thd & dists ~= 0);
        
        if ~isempty(ind)
            map.pw_pos_matrix_plus(indx(ind), indy(ind)) = 0;
            [indx, indy] = find(map.pw_pos_matrix_plus > 0);
            ii = 1; %% set ii back
            
        else
            ii = ii + 1;
        end
    end
    
    %%% Clean up negatively charged pinwheels
    [indx, indy] = find(map.pw_pos_matrix_minus > 0);
    ii = 1;
    
    while ii < length(indx)
        dists = sqrt((indx(ii) - indx).^2 + (indy(ii) - indy).^2);        
        ind = find(dists < thd & dists ~= 0);
        if ~isempty(ind)
            map.pw_pos_matrix_minus(indx(ind), indy(ind)) = 0;
            [indx, indy] = find(map.pw_pos_matrix_minus > 0);
            ii = 1; %% set ii back
        else
            ii = ii + 1;
        end
    end
    
    [PWposY_plus, PWposX_plus] = find(map.pw_pos_matrix_plus > 0);
    [PWposY_minus, PWposX_minus]  = find(map.pw_pos_matrix_minus > 0);
    % NOTE: x and y indices are reversed in the find command so Y here is
    % first and x is second!!! 
        

    if ~isdeployed && do_plotting
        load colormap_min;
        
        figure()
        imagesc(angle(map.z_filtered));
        colormap(colormap_min);
        set(gca, 'ydir','normal');
        hold on;
        plot(PWposX_plus,PWposY_plus,'o','color','w','MarkerSize',7, 'linewidth', 2);
        plot(PWposX_minus, PWposY_minus,'o','color','k','MarkerSize',7, 'linewidth', 2);
        set(gca, 'ydir', 'normal');
        axis image;
        hold off;
        
        figure;
        imagesc(map.automated_pw_density);    
        set(gca, 'ydir','normal');
        axis image;
        colormap(jet);
        colorbar;
        hold on;
        plot(PWposX_plus,PWposY_plus,'o','color','w','MarkerSize',10);
        plot(PWposX_minus, PWposY_minus,'o','color','w','MarkerSize',10);
        hold off;
        title('Automated pw dens estimates and pw positions');
    end
    
    
    %%%% Store in map structure
    map.PWxList = [PWposX_plus; PWposX_minus];
    map.PWyList = [PWposY_plus ; PWposY_minus];
    map.signlist = [ones(size(PWposX_plus)) ; -ones(size(PWposX_minus))];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Pinwheel density with pw pos estimates is: ' num2str(length(map.PWxList)/(sum(map.new_roi(:))/map.average_w^2))]);
    disp(['Pinwheel density with plateu fitting  is: ' num2str(lpwd)]);

    
end %% END OF FUNCTION ESTIMATE_LOCAL_PW_DENSITY.M
