function [map, absolute_scale, measure,llp,hhp,llp_cutoffs] = filter_orientation_map(map,shrinkage, measure,absolute_scale,llp,hhp,llp_cutoffs,beta, do_plotting)
%
% FUNCTION CALLED BY ANALYZE_SINGLE_MAP.M
%
% USAGE : 
% 
% [map, absolute_scale, measure,llp,hhp,llp_cutoffs] = ...
% filter_orientation_map(map,shrinkage, measure,absolute_scale,llp,hhp,llp_cutoffs,beta, do_plotting)
%
% INPUT PARAMETERS:
% map           ... data structure containing all information (see help ANALYZE_SINGLE_MAP.M for more info
% shrinkage     ... choose 1 of no shrinkage desired, choose 2 if
% measure       ... this parameter determines how many pixels correspond to
%                   one millimeter, depending on shrinkage, this might
%                   change and the new return value for measure should be
%                   used subsequently for processing
%
% absolute scale... this parameter determines how millimeters correspond to one pixel, depending on shrinkage, this might
%                   change and the new return value for absolute should be
%                   used subsequently for processing
%  llp         ... lowpass cutoff frequency in pixels, depending on shrinkage, this might
%                   change and the new return value for llp should be
%                   used subsequently for processing
%  hhp         ... highpass cutoff frequency in pixels, depending on shrinkage, this might
%                   change and the new return value for hhp should be
%                   used subsequently for processing
%  beta         ... "Temperature" for the fermi-filter if used,
%                   beta is measured in units of k_hp or k_lp (see Kaschube et al., Science 2010)
%                   choose beta = 0 for gaussian filtering (not recommended)
%
%
% do_plotting   ... flag, either zero or one, chose 1 for visual output of
%                   filter results
%
%  Copyright (C) 2014, Max-Planck-Institut for Dynamics and Self-organization, The  
%  Nonlinear Dynamics Group. This software may be used, copied, or 
%  redistributed as long as it is not sold and this copyright notice is 
%  reproduced on each copy made. This routine is provided as is without 
%  any express or implied warranties whatsoever.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %%%%%%%%%%%%%%%%%%%%% Shrinkage is done, only if desired 
    map.shrinkage = shrinkage;
    if shrinkage > 1
        disp(['Image is shrinked by a factor of ' num2str(shrinkage)]);
        %%%% Shrink image (faster processing and to avoid memory overflow!)
        data = shrink(map.z_raw,shrinkage);
        roi_matrix = shrink(map.roi_raw,shrinkage);
        roi_matrix(roi_matrix < 1) = 0;
        measure = measure/shrinkage;
        absolute_scale = absolute_scale*shrinkage;

        hhp = hhp/shrinkage;
        llp = llp/shrinkage;
        llp_cutoffs = llp_cutoffs./shrinkage;
        %%%%%%%%%%%%%%%% map.z_shrinked contains raw data without any filtering
        map.z_shrunk = data;
        map.roi_shrunk = roi_matrix;
        
    else %%% No shrinking
        
        %%%%%%%%%%%%%%%% map.z_shrinked contains raw data without any filtering
        map.z_shrunk = map.z_raw;
        map.roi_shrunk = map.roi_raw;

    end
    
    
    
    if ~isdeployed && do_plotting
       scrsz = get(0,'ScreenSize');

        figure(1);
        set(gcf,'Position',[1 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3]); %[left, bottom, width, height]:
        subplot(2,2,1);
        imagesc(real(map.z_shrunk))
        set(gca, 'YDir', 'normal');
        colormap(gray);
        axis image;
        title('Unprocessed real part');
        subplot(2,2,2);
        imagesc(imag(map.z_shrunk))
        set(gca, 'YDir', 'normal');
        colormap(gray);
        axis image;
        title('Unprocessed imaginary part');

    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Now do Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    data_tmp = map.z_shrunk; 
    
    %%% Set everything outside of the ROI to zero
    data_tmp(map.roi_shrunk == 0) = 0;


    disp('Doing Highpass Fermi Filtering...')
    % Highpass
    [data_tmp]= get_filtered_for_complex_fields(data_tmp,map.roi_shrunk,'highpass',hhp,beta);
    %%% Set everything outside of the ROI to zero
    data_tmp(map.roi_shrunk == 0) = 0;
    disp('...done');

    % Lowpass
    disp('Doing Lowpass Fermi Filtering...')
    [data_tmp]= get_filtered_for_complex_fields(data_tmp,map.roi_shrunk, 'lowpass',llp,beta);
    %%% Set everything outside of the ROI to zero
    data_tmp(map.roi_shrunk == 0) = 0;
    
    %%%% Set variance to one
    data_tmp(map.roi_shrunk == 1) = data_tmp(map.roi_shrunk == 1)./std(data_tmp(map.roi_shrunk == 1));
    disp('...done');
        
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    map.z_filtered = data_tmp;   

    
    %%%% DO PLOTTING IF DESIRED
    
    if ~isdeployed && do_plotting
        load colormap_min;
        figure(1);
        subplot(2,2,3);
        imagesc(real(map.z_filtered));
        colormap(gray);
        set(gca, 'YDir', 'normal');
        colormap(gray);
        axis image;
        title('Real part of bandpass filtered data');    
        subplot(2,2,4);
        imagesc(imag(map.z_filtered));
        colormap(gray);
        set(gca, 'YDir', 'normal');
        colormap(gray);
        axis image;
        title('Imaginary part of bandpass filtered data');    

        figure;
        set(gcf,'Position',[scrsz(3)/3 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3]); %[left, bottom, width, height]:
        imagesc(angle(map.z_filtered));
        axis image;
        colormap(colormap_min);
        set(gca, 'ydir', 'normal');
        title('Angle map of bandpass filtered data');
        drawnow
           
    end
    
    
    
end %%%% END OF FUNCTION FILTER_ORIENTATION_MAP.M