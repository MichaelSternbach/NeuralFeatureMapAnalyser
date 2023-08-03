function [pinwheel_stats,pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track)
%{
 [pinwheel_stats,pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track)

 Code that calculates the corresponding pinwheels of the map in multiple
 bootstrap samples. As the program runs a progress report is displayed. To
 prevent this from happening check do_report after this definition.

 INPUT
 data_obj: handle to class that does the data handling, eg. extracting
       bootstrap samples, reading ROI, filtering, etc.
       => see class: data_handle
 tracker_obj: handle to class that does the tracking between input maps
       => see class: pinwheel_tracker
 Set the properties in the corresponding classes!
 simple_track: compare only MEAN sample with the rest, not across sub-samples 

 OUPUT:
 pinwheel_stats: structure with fields with tables corresponding to
       matching pinwheels in the different bootstrap samples
       - x = position in X 
       - y = position in Y
       - sign = sign
       - label = number assigned to pinwheel when tracker_obj pinwheel
              finder class is used to find pinwheels
       - probability = how often each pinwheel is found between samples
       Fields filled with NAN mean there is no matching pinwheel.
 pinwheel_spurious: same structure, but contains a list of pinwheels that
                   are not included in pinwheel_output. Useful when doing
                   simple_track

 chepe@nld.ds.mpg.de   3/Nov/2015
%}
%%

% check input
if nargin==2
    simple_track = false;
end

% do a report of the progress?
do_report = true;

% read the number of samples to be used
boot_samples = size(data_obj.samples_array,3);

%% Create positions array using the first sample as base
% Make cell structures where the matching pinwheels between comparisons are
% stored. The first cell is the information of the pinwheels for the sample.
% The second cell has the comparison between elements. 1st element = sample1 vs rest
% 2nd element = sample2 vs rest, etc.
% Starting from the 2nd, only the pinwheels not found in the 1st are stored
% and compared.
pinwheel_tables = cell(boot_samples,1);
tracking_tables = cell(boot_samples,1);

%% First round: compare 1st sample vs the rest

% get pinwheels in 1st sample and fill tracking tables
base = 1;
z_base = data_obj.filter_map(data_obj.read_map(base));
pinwheel_tables{base} = tracker_obj.find_pinwheels(z_base,data_obj.ROI,base);
tracking_tables{base} = zeros(length(pinwheel_tables{base}.number),boot_samples);
tracking_tables{base}(:,base) = 1:length(pinwheel_tables{base}.number);

% start comparing 1st sample against the rest (labeled as input)
if do_report
disp(['Doing comparisons with base = 1 num_pinwheels = ',num2str(length(pinwheel_tables{1}.number))])
end

for input=2:boot_samples
    z_input = data_obj.filter_map(data_obj.read_map(input));
    
    % get pinwheels in input sample and fill tracking tables
    pinwheel_tables{input}=tracker_obj.find_pinwheels(z_input,data_obj.ROI,base);
    
    if do_report
    time = clock;
    disp(['    time = ',num2str(time(4)),':',num2str(time(5),'%02d'),' input = ',num2str(input),' num_pinwheels = ',num2str(length(pinwheel_tables{input}.number))])
    end
    
    % track the interpolation between the two maps
    tracking = tracker_obj.interpolate(z_base,z_input,data_obj.ROI);
    
    % fill the tracking tables
    if ~isempty(tracking.ini)
        % --> corresponding pinwheels between 1st sample and input
        tracking_tables{1}(:,input)=tracking.ini(:,2);
    end
    if ~isempty(tracking.end)
        % --> pinwheels in input that are not part of 1st sample: add to
        % new index in the tracking_table
        tracking_tables{input}=zeros(size(tracking.end,1),boot_samples);
        tracking_tables{input}(:,[1 input])=tracking.end;
    end
end

%% If simple track done, gather results and return
if simple_track
    
    % Make empty array
    % [pinwheel, sample, [number,x,y,sign]]
    pinwheel_stats = struct(...
        'x',NaN*zeros(size(pinwheel_tables{1}.x,1),boot_samples),...
        'y',NaN*zeros(size(pinwheel_tables{1}.x,1),boot_samples),...
        'sign',NaN*zeros(size(pinwheel_tables{1}.x,1),boot_samples),...
        'label',NaN*zeros(size(pinwheel_tables{1}.x,1),boot_samples),...
        'probability',zeros(size(pinwheel_tables{1}.x,1),1));
       
    for sampleNum=1:boot_samples
        % find indices in pinwheel data 
        pw_ind = tracking_tables{1}(:,sampleNum);
        % fill the table
        pinwheel_stats.label(:,sampleNum) = pw_ind;
        pinwheel_stats.x(pw_ind>0,sampleNum) = pinwheel_tables{sampleNum}.x(pw_ind(pw_ind>0));
        pinwheel_stats.y(pw_ind>0,sampleNum) = pinwheel_tables{sampleNum}.y(pw_ind(pw_ind>0));
        pinwheel_stats.sign(pw_ind>0,sampleNum) = pinwheel_tables{sampleNum}.sign(pw_ind(pw_ind>0));  
    end
    pinwheel_stats.probability = sum(pinwheel_stats.label>0,2)/boot_samples;
      
    % make extra output of pinwheels that are not matched from the other tables
    if nargout == 2
        pinwheel_spurious = struct('x',[],'y',[],'sign',[],'number',[],'label',[]);
        
        for ii=2:length(tracking_tables)
            if isempty(tracking_tables{ii})
                continue
            end
            
            ID_table = tracking_tables{ii}(:,ii);
            ID_pinwheel = pinwheel_tables{ii}.number;
            
            [~,test] = intersect(ID_pinwheel,ID_table);
            
            pinwheel_spurious.x = [pinwheel_spurious.x;pinwheel_tables{ii}.x(test)];
            pinwheel_spurious.y = [pinwheel_spurious.y;pinwheel_tables{ii}.y(test)];
            pinwheel_spurious.sign = [pinwheel_spurious.sign;pinwheel_tables{ii}.sign(test)];
            pinwheel_spurious.number = [pinwheel_spurious.number;pinwheel_tables{ii}.number(test)];
            tmp = pinwheel_tables{ii}.number(test);
            pinwheel_spurious.label = [pinwheel_spurious.label;[ii*ones(size(tmp)) tmp]];
            
        end
    end

    return
end
%% 2nd-Nth round: compare 2nd vs rest, 3rd vs rest, etc
% NOTE: If a matching pinwheel between two tracking_tables are found, then
% that pinwheel in the input table is removed. This process is continued
% until all comparisons are done or the tracking tables are empty

% Add first table to final combined table
combined_table = tracking_tables{1};

% Use different maps as base for comparison
ranked_samples = 2:boot_samples;
while ~isempty(ranked_samples)
    
    % Sort next base map order to make less comparissons
    % Rank samples depending on similarity with first sample to sort the order
    % in which the pinwheels are matched. Samples where there is nothing to
    % match are removed
    % lost = sum(tracking_tables{1}(ranked_samples,:)<1,1);
    new = zeros(size(ranked_samples));
    for ind=1:length(ranked_samples)
        new(ind)=size(tracking_tables{ranked_samples(ind)},1);
    end
    new(new==0)=NaN;
    [new,ind]=sort(new,'descend');
    ranked_samples = ranked_samples(ind);
    ranked_samples(isnan(new))=[];
    
    % check if all samples are complete
    if isempty(ranked_samples)
        break
    end
    
    % get base sample number and map
    base = ranked_samples(1);
    z_base = data_obj.filter_map(data_obj.read_map(base));

    % Use different maps as input for comparison
    if do_report
    total = 0;
    for ii=2:length(ranked_samples)
        total = total + size(tracking_tables{ranked_samples(ii)},1);
    end
    disp(['Doing comparisons with base = ',num2str(base),', missing pinwheels = ',num2str(size(tracking_tables{base},1)),', total to match = ',num2str(total)])
    end
    
    % Match pinwheels between samples
    for input_ind=2:length(ranked_samples)
        if do_report
        time = clock;
        disp(['    time = ',num2str(time(4)),':',num2str(time(5),'%02d'),', input = ',num2str(ranked_samples(input_ind)),', missing pinwheels = ',num2str(size(tracking_tables{ranked_samples(input_ind)},1))])
        end
        
        % interplate between the base map and the input map
        input = ranked_samples(input_ind);
        z_input = data_obj.filter_map(data_obj.read_map(input));
        
        % track the interpolation between the two maps
        tracking = tracker_obj.interpolate(z_base,z_input,data_obj.ROI);
        
        % fill tracking tables with results and SORT to use intersect
        base_list = sortrows(tracking.ini,1);
        input_list = sortrows([tracking.ini(tracking.ini(:,2)>0,:);tracking.end],2);
        
        [~,ind_result,ind_table]=intersect(base_list(:,1),tracking_tables{base}(:,base));
        tracking_tables{base}(ind_table,input) = base_list(ind_result,2);
        
        [~,ind_result,ind_table]=intersect(input_list(:,2),tracking_tables{input}(:,input));
        tracking_tables{input}(ind_table,base) = input_list(ind_result,1);
        
        % Delete found matching pinwheels from input tracking_table
        ind_delete = ismember(tracking_tables{input}(:,input),tracking_tables{base}(:,input));
        tracking_tables{input}(ind_delete,:)=[];
    end
    
    % before adding, set matched pinwheels in base that are already part of
    % combined_table to -Inf
    for col_num=1:boot_samples
        %[~,idx_row] = intersect(tracking_tables{base}(:,col_num),combined_table(  combined_table(:,col_num)>0  ,col_num));
        idx_row = ismember(tracking_tables{base}(:,col_num),combined_table(  combined_table(:,col_num)>0  ,col_num));
        tracking_tables{base}(idx_row,col_num) = -Inf;
    end
    
    % add filled base table to combined_table
    combined_table = [combined_table;tracking_tables{base}];
    
    % remove current base index
    ranked_samples(1)=[];
    
end

%% Sanity check: are pinwheels repeating/missing??

for sample2test=1:boot_samples
    num_pw = length(pinwheel_tables{sample2test}.number);
    test_input=sort(combined_table(combined_table(:,sample2test)>0,sample2test));
    if sum(diff(test_input)==0)
        error(['Problem when combining tracking tables in sample ',num2str(sample2test),': Repeated pinwheels'])
    end
    if length(unique(test_input))~=num_pw
        error(['Problem when combining tracking tables in sample ',num2str(sample2test),': Missing pinwheels'])
    end
end

%% Prepare output

% [pinwheel, sample, [number,x,y,sign]]
pinwheel_stats = struct(...
    'x',NaN*zeros(size(combined_table,1),size(combined_table,2)),...
    'y',NaN*zeros(size(combined_table,1),size(combined_table,2)),...
    'sign',NaN*zeros(size(combined_table,1),size(combined_table,2)),...
    'label',NaN*zeros(size(combined_table,1),size(combined_table,2)),...
    'probability',zeros(size(combined_table,1),1));

for sampleNum=1:boot_samples
    % find indices in pinwheel data
    pw_extracted = pinwheel_tables{sampleNum};
    ind_table =find(combined_table(:,sampleNum)>0);
    ind_pw =combined_table(combined_table(:,sampleNum)>0,sampleNum);
    
    % fill data
    pinwheel_stats.x(ind_table,sampleNum)=pw_extracted.x(ind_pw);
    pinwheel_stats.y(ind_table,sampleNum)=pw_extracted.y(ind_pw);
    pinwheel_stats.sign(ind_table,sampleNum)=pw_extracted.sign(ind_pw); % min(pinwheel_output.sign(ind_table),pw_extracted.sign(ind_pw)); % min() to remove NaN
    pinwheel_stats.label(:,sampleNum)=combined_table(:,sampleNum);
    pinwheel_stats.probability(ind_table)=pinwheel_stats.probability(ind_table)+~isnan(pw_extracted.x(ind_pw))/boot_samples;
end

% since all the tables have been joined, there are no remaining pinwheels
pinwheel_spurious = [];


end
