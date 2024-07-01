function map =makeMapRad(data,stimuli_cond)
%     map = zeros(size(data,1:2));
%     for ii_condition = 1:size(stimuli_cond)
%         if isnan(stimuli_cond(ii_condition))
%             continue
%         end
%         tmp_img = mean(data(:,:,ii_condition,:),4);
%         map = map + exp(1i*real(stimuli_cond(ii_condition)))*tmp_img;
%     end

    map = sum(data.*exp(1i*real(reshape(stimuli_cond,[1 1 size(stimuli_cond,1)]))),3);
end
