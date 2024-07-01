function data_obj_rand = getRandomisedDataHandle(data_obj)
    data = randomizeData(data_obj.data);
    data_obj_rand =  data_handle_corrected(data_obj.info,data,data_obj.ROI);
end

function data_rand = randomizeData(data)
    data_rand = zeros(size(data));
    rand_order = shuffle(1:(size(data,3)*size(data,4)));
    for ii_stim = 1:size(data,3)
        for ii_trial = 1:size(data,4)
            rand_num=rand_order((ii_stim-1)*size(data,4)+ii_trial);
%             ii_trial_rand = mod(rand_num,size(data,4));
%             ii_stim_rand = (rand_num-ii_trial_rand)/size(data,4);
            [row,col] = ind2sub(size(data,[3 4]),rand_num);
            data_rand(:,:,ii_stim,ii_trial) = data(:,:,row,col);
        end
    end

end

%% from https://de.mathworks.com/matlabcentral/answers/20789-shuffle-numbers-in-a-vector
 function v=shuffle(v) 
     v=v(randperm(length(v)));
 end