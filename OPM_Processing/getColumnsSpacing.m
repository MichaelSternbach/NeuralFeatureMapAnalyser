function [average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient,CI_average_spacing_mm,CI_local_spacing_mm] = getColumnsSpacing(data_obj,DataFolder,smallest_w_mm,largest_w_mm,w_step_mm,getCI,FilterMap,alpha,direction_data)
    if nargin < 6
        getCI = false;
    end
    if nargin < 7
        FilterMap = false;
    end
    if nargin < 8
        alpha = 0.05;
    end
    if nargin < 9
        direction_data = false;
    end

    %% get mean spacing
    if direction_data
        dir_str = '_DirectionMap';
    else
        dir_str = '';
    end

    if FilterMap
        SpacingFile = [DataFolder 'MapSpacingFiltered_' data_obj.info.ID dir_str '.mat'];
    else
        SpacingFile = [DataFolder 'MapSpacing_' data_obj.info.ID dir_str '.mat'];
    end
    if isfile(SpacingFile)
        disp(['load spacing data from: ' SpacingFile])
        load(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI','WavletCoefficient')
    else
        if FilterMap
            z = data_obj.filter_map(data_obj.read_map(1,false,direction_data));
        else
            z = data_obj.read_map(1,false,direction_data);
        end
        [average_spacing_mm,local_spacing_mm,newROI,WavletCoefficient] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
        save(SpacingFile,'average_spacing_mm','local_spacing_mm','newROI','WavletCoefficient')
    end
    
    %% get CIs spacing for local and mean column spacing
    % calculate column spacing for all bootstrap samples
    % and all Jackknife samples
    % based on their values calculate confidence intervals via Bca Method
    % described in Efron - Computer Age Statistical Inference
    if nargin < 8
        alpha = 0.05;   
    end

    if getCI == 1
        if FilterMap
            CISpacingFile = [DataFolder 'CI_MapSpacingFiltered_' data_obj.info.ID dir_str '.mat'];
        else
            CISpacingFile = [DataFolder 'CI_MapSpacing_' data_obj.info.ID dir_str '.mat'];
        end
        if isfile(CISpacingFile)
            
            FileData = load(CISpacingFile,'alpha');
            if FileData.alpha == alpha
                disp('load ColumnSpacingCIs')
                load(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm','average_spacings_mm','local_spacings_mm','newROIs')
                disp('data loaded !')
            else
                disp('Filedata has different alpha!')
                
                disp('load BS and JS')
                load(CISpacingFile,'average_spacings_mm','local_spacings_mm','jackstat_average_spacing_mm','local_spacingsJS_mm','newROIsBS','newROIsJS')
                

                bootstat_local_spacings_mm = convert2vector(local_spacings_mm,data_obj,2);
                jackstat_local_spacing_mm =  convert2vector(local_spacingsJS_mm,data_obj,1);

                disp('calc CS CIs fro BS and JS')
                CI_average_spacing_mm = bootstrap_ci(average_spacings_mm,average_spacing_mm,jackstat_average_spacing_mm,alpha);
                CI_local_spacing_mmVector = bootstrap_ci(bootstat_local_spacings_mm,data_obj.array2vector(local_spacing_mm),jackstat_local_spacing_mm,alpha);           
                CI_local_spacing_mm = zeros([size(data_obj.ROI) 2]);
                CI_local_spacing_mm(:,:,1) = data_obj.vector2array(CI_local_spacing_mmVector(:,1));
                CI_local_spacing_mm(:,:,2) = data_obj.vector2array(CI_local_spacing_mmVector(:,2));
            end

            
        else
            num_boot_samples = size(data_obj.samples_array,3);
            %% get spacing of bootstraped map 
            average_spacings_mm = zeros(1,size(data_obj.samples_array,3));
            bootstat_local_spacings_mm = zeros(sum(data_obj.ROI(:)),num_boot_samples);
            local_spacings_mm = cell(1,num_boot_samples);
            newROIsBS = cell(1,num_boot_samples);
            parfor ii = 2:num_boot_samples
                disp(['BS' num2str(ii)])
                if FilterMap
                    z = data_obj.filter_map(data_obj.read_map(ii,false,direction_data));
                else
                    z = data_obj.read_map(ii,false,direction_data);
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
            local_spacingsJS_mm = cell(1,data_obj.data_parameters.num_blocks);
            newROIsJS = cell(1,data_obj.data_parameters.num_blocks);
            parfor ii=1:data_obj.data_parameters.num_blocks
                disp(['JS' num2str(ii)])
                if FilterMap
                    z = data_obj.filter_map(data_obj.read_map(ii,false,direction_data));
                else
                    z = data_obj.read_map(ii,false,direction_data);
                end
                [average_spacing_mm_js,local_spacing_mm_js,newROI_js] = get_column_spacingManuel(z,data_obj.ROI,data_obj.info.pix_per_mm,smallest_w_mm,largest_w_mm,w_step_mm);
                
                jackstat_average_spacing_mm(ii) = average_spacing_mm_js;
                jackstat_local_spacing_mm(:,ii) = data_obj.array2vector(local_spacing_mm_js);
                local_spacingsJS_mm{ii} = local_spacing_mm_js;
                newROIsJS{ii} = newROI_js; 
            end
            data_obj.set_samples_array(samples_array);

            
            %% calculate confidence intervals
            CI_average_spacing_mm = bootstrap_ci(average_spacings_mm,average_spacing_mm,jackstat_average_spacing_mm,alpha);
            
            CI_local_spacing_mmVector = bootstrap_ci(bootstat_local_spacings_mm,data_obj.array2vector(local_spacing_mm),jackstat_local_spacing_mm,alpha);           
            CI_local_spacing_mm = zeros([size(data_obj.ROI) 2]);
            CI_local_spacing_mm(:,:,1) = data_obj.vector2array(CI_local_spacing_mmVector(:,1));
            CI_local_spacing_mm(:,:,2) = data_obj.vector2array(CI_local_spacing_mmVector(:,2));
            
            save(CISpacingFile,'CI_average_spacing_mm','CI_local_spacing_mm','average_spacings_mm','local_spacings_mm','jackstat_average_spacing_mm','local_spacingsJS_mm','newROIsBS','newROIsJS','alpha')
        end
    else
        CI_average_spacing_mm = [];
        CI_local_spacing_mm = [];
    end
end


function VectorArray = convert2vector(ArrayCell,data_obj,first)
    for ii = first:length(ArrayCell)
        vector = data_obj.array2vector(ArrayCell{ii});
        if ii == first
            VectorArray = zeros([length(vector) length(ArrayCell)]);
        end
        VectorArray(:,ii) = vector;
    end
end