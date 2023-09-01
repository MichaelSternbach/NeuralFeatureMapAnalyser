function CalcHyperColumnSizes(animal,experiment_numList,smallest_w,w_step,folder)
    % smallest_w = 0.5;
    % w_step = 0.05;
    % largest_w = 1.5;
    % 
    % animal = 'Mouse Lemur';
    % experiment_numList = 1:5;
    % folder = ''

    for ii = 1:size(experiment_numList,2)
        experiment_num = experiment_numList(ii);
        [data_info,~,data_obj,~,~] = getAnimalData(animal,experiment_num,1);

        z_input = data_obj.filter_map(data_obj.read_map(1));

        SpacingFile = [folder 'MapSpacing_' animal num2str(experiment_num) '.mat'];
%         if isfile(SpacingFile)
%             load(SpacingFile,'average_spacing_mm','local_spacing_mm')
%         else
%             [average_spacing_mm,local_spacing_mm] = get_column_spacingManuel(z_input,data_obj.ROI,data_info.pix_per_mm,smallest_w,largest_w,w_step);
%             save(SpacingFile,'average_spacing_mm','local_spacing_mm')
%         end
        if ~isfile(SpacingFile)
            [average_spacing_mm,local_spacing_mm] = get_column_spacingManuel(z_input,data_obj.ROI,data_info.pix_per_mm,smallest_w,largest_w,w_step);
            save(SpacingFile,'average_spacing_mm','local_spacing_mm')
        end



    end

end


