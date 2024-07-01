function [average_spacing_mm,local_spacing_mm,newROI,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap)
    if nargin < 6
        getCI = false;
    end
    if nargin < 7
        FilterMap = true;
    end

    %% get mean spacing
    SpacingFile = [DataFolder 'MapSpacing_' data_obj.info.ID '.mat'];
    if isfile(SpacingFile)
        load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI')
    else
        if FilterMap
            z = data_obj.filter_map(data_obj.read_map());
        else
            z = data_obj.read_map();
        end
        [average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
        save(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI','WavletCoefficient')
    end
    
    if nargin < 6
        getCI = false;
    end
    %% get CIs
    alpha = 0.05;
    if getCI == 1
        CISpacingFile = [DataFolder 'CI_MapSpacing_' data_obj.info.ID '.mat'];
        if isfile(CISpacingFile)
            load(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm','average_spacings_mm','local_spacings_mm','newROIs')
        else
            num_boot_samples = size(data_obj.samples_array,3);
            %% get spacing of bootstraped map 
            average_spacings_mm = zeros(1,size(data_obj.samples_array,3));
            bootstat_local_spacings_mm = zeros(sum(data_obj.ROI(:)),num_boot_samples);
            for ii = 2:num_boot_samples
                if FilterMap
                    z = data_obj.filter_map(data_obj.read_map(ii));
                else
                    z = data_obj.read_map(ii);
                end
                [average_spacing_mm_bs,local_spacing_mm_bs,newROI_bs] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
                average_spacings_mm(ii) = average_spacing_mm_bs;
                local_spacings_mm{ii} = local_spacing_mm_bs;
                newROIsBS{ii} = newROI_bs; 
                
                bootstat_local_spacings_mm(:,ii) = data_obj.array2vector(local_spacing_mm_bs);
                
            end
            
            %% get spacing of jackknife samples
            samples_array = data_obj.samples_array;
            data_obj.prepare_jackknife_samples;
            jackstat_average_spacing_mm = zeros(1,data_obj.data_parameters.num_blocks);
            jackstat_local_spacing_mm = zeros(sum(data_obj.ROI(:)),data_obj.data_parameters.num_blocks);
            for ii=1:data_obj.data_parameters.num_blocks
                if FilterMap
                    z = data_obj.filter_map(data_obj.read_map(ii));
                else
                    z = data_obj.read_map(ii);
                end
                [average_spacing_mm_js,local_spacing_mm_js,newROI_js] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
                
                jackstat_average_spacing_mm(ii) = average_spacing_mm_js;
                jackstat_local_spacing_mm(:,ii) = data_obj.array2vector(local_spacing_mm_js);
                local_spacingsJS_mm{ii} = local_spacing_mm_js;
                newROIsJS{ii} = newROI_js; 
            end
            data_obj.set_samples_array(samples_array);
            
            %% 
            CI_average_spacing_mm = bootstrap_ci(average_spacings_mm,average_spacing_mm,jackstat_average_spacing_mm,alpha);
            
            CI_local_spacing_mmVector = bootstrap_ci(bootstat_local_spacings_mm,data_obj.array2vector(local_spacing_mm),jackstat_local_spacing_mm,alpha);           
            CI_local_spacing_mm = zeros([size(data_obj.ROI) 2]);
            CI_local_spacing_mm(:,:,1) = data_obj.vector2array(CI_local_spacing_mmVector(:,1));
            CI_local_spacing_mm(:,:,2) = data_obj.vector2array(CI_local_spacing_mmVector(:,2));
            
            save(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm','average_spacings_mm','local_spacings_mm','jackstat_average_spacing_mm','local_spacingsJS_mm','newROIsBS','newROIsJS','alpha')
        end
    end
end