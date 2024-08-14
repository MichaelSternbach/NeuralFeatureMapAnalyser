function data_rand = randomizeData(data)
    data_rand = zeros(size(data));
    rand_order = shuffle(1:(size(data,3)*size(data,4)));
    stim_list = zeros(1,size(data,3)*size(data,4));
    stim_list_rand = zeros(1,size(data,3)*size(data,4));
%     for ii_stim = 1:size(data,3)
%         for ii_trial = 1:size(data,4)
%             rand_num=rand_order((ii_stim-1)*size(data,4)+ii_trial);
% %             ii_trial_rand = mod(rand_num,size(data,4));
% %             ii_stim_rand = (rand_num-ii_trial_rand)/size(data,4);
%             [row,col] = ind2sub(size(data,[3 4]),rand_num);
%             data_rand(:,:,ii_stim,ii_trial) = data(:,:,row,col);
% 
%             stim_list((ii_stim-1)*size(data,4)+ii_trial) = ii_stim;
%             stim_list_rand((ii_stim-1)*size(data,4)+ii_trial) = row;
%         end
%     end
    for ii = 1:(size(data,3)*size(data,4))
        rand_num=rand_order(ii);
        [row,col] = ind2sub(size(data,[3 4]),ii);
        [row_rand,col_rand] = ind2sub(size(data,[3 4]),rand_num);
        data_rand(:,:,row,col) = data(:,:,row_rand,col_rand);

        stim_list(ii) = row;
        stim_list_rand(ii) = row_rand;
    end
    % index = 1:(size(data,3)*size(data,4));
    % figure;plot(index,rand_order,'*')
    % figure;histogram2(stim_list,stim_list_rand)
    % figure;plot(stim_list,stim_list_rand,'*')
end

%% from https://de.mathworks.com/matlabcentral/answers/20789-shuffle-numbers-in-a-vector
 function v=shuffle(v) 
     v=v(randperm(length(v)));
 end