function [pinwheel_stats,pinwheel_spurious] = trackPws3D(z_filtered,ROI,tracker_obj)
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


% read the number of samples to be used
depth = size(z_filtered,1);

% do a report of the progress?
do_report = true;



%% Create positions array using the first sample as base
% Make cell structures where the matching pinwheels between comparisons are
% stored. The first cell is the information of the pinwheels for the sample.
% The second cell has the comparison between elements. 1st element = sample1 vs rest
% 2nd element = sample2 vs rest, etc.
% Starting from the 2nd, only the pinwheels not found in the 1st are stored
% and compared.
pinwheel_tables = cell(depth+1,1);
tracking_tables = cell(depth+1,1);

%% First round: compare 1st sample vs the rest

% % get pinwheels in 1st sample and fill tracking tables
% 
% pinwheel_tables{base} = tracker_obj.find_pinwheels(z_base,data_obj.ROI,base);
% tracking_tables{base} = zeros(length(pinwheel_tables{base}.number),depth);
% tracking_tables{base}(:,base) = 1:length(pinwheel_tables{base}.number);

% start comparing 1st sample against (labeled as input)
base = 1;

% z1 = reshape(z_filtered(base,:,:),size(z_filtered,[2 3]));
% ROI2D = reshape(ROI(base,:,:),size(ROI,[2 3]));
% 
% pinwheel_tables{base} = tracker_obj.find_pinwheels(z1,ROI2D,base);
% tracking_tables{base} = zeros(length(pinwheel_tables{base}.number),depth);
% tracking_tables{base}(:,base) = 1:length(pinwheel_tables{base}.number);

ROI2D = reshape(sum(ROI,1),size(ROI,[2 3]));%reshape(ROI(jj,:,:),size(ROI,[2 3])).*reshape(ROI(jj-1,:,:),size(ROI,[2 3]));
% ROI2D = reshape(ROI(1,:,:),size(ROI,[2 3]));
ROI2D = (ROI2D>0);

samples = (base+1):depth;

for jj= samples
    
%     ROI2D = reshape(ROI(jj,:,:),size(ROI,[2 3])).*reshape(ROI(jj-1,:,:),size(ROI,[2 3]));
%     ROI2D = (ROI2D==1);
    if jj == base+1
        %z1 = reshape(z_filtered(jj-1,:,:),size(z_filtered,[2 3]));
        z1 = reshape(z_filtered(base,:,:),size(z_filtered,[2 3]));

        pinwheel_tables{base} = tracker_obj.find_pinwheels(z1,ROI2D,base);
        tracking_tables{base} = zeros(length(pinwheel_tables{base}.number),depth);
        tracking_tables{base}(:,base) = 1:length(pinwheel_tables{base}.number);
        if do_report
            disp(['Finished calculations for base =' num2str(base) 'num_pinwheels = ',num2str(length(pinwheel_tables{base}.number))])
        end
    else
        z1 = z2;
    end
    z2 = reshape(z_filtered(jj,:,:),size(z_filtered,[2 3]));
    
    
    % get pinwheels in input sample and fill tracking tables
    pinwheel_tables{jj}=tracker_obj.find_pinwheels(z2,ROI2D,jj);
    
    if do_report
    time = clock;
    disp(['    time = ',num2str(time(4)),':',num2str(time(5),'%02d'),' input = ',num2str(jj),' num_pinwheels = ',num2str(length(pinwheel_tables{jj}.number))])
    end
    
    % track the interpolation between the two maps
    [tracking,tracking_table,tracking_structure] = tracker_obj.interpolate(z1,z2,ROI2D);
    tracking_tables{jj}=tracking;
    
    plotPinwheels(z1,ROI2D,pinwheel_tables{jj-1},tracking.ini(:,1))
    plotPinwheels(z2,ROI2D,pinwheel_tables{jj},tracking.ini(:,2))
    close all

%     % fill the tracking tables
%     if ~isempty(tracking.ini)
%         % --> corresponding pinwheels between 1st sample and input
%         tracking_tables{base}(:,jj)=tracking.ini(:,2);
%     end
%     if ~isempty(tracking.end)
%         % --> pinwheels in input that are not part of 1st sample: add to
%         % new index in the tracking_table
%         tracking_tables{jj}=zeros(size(tracking.end,1),depth);
%         tracking_tables{jj}(:,[1 jj])=tracking.end;
%     end
end

%% If simple track done, gather results and return


% Make empty array
% [pinwheel, sample, [number,x,y,sign]]
Size_Pw = size(pinwheel_tables{base}.x,1)*depth;
pinwheel_stats = struct(...
    'x',NaN*zeros(Size_Pw,depth),...
    'y',NaN*zeros(Size_Pw,depth),...
    'sign',NaN*zeros(Size_Pw,depth),...
    'label',zeros(Size_Pw,depth));
   
for sampleNum=1:depth
    % find indices in pinwheel data 
%     pw_ind = tracking_tables{sampleNum}(:,sampleNum);
    if sampleNum ==1
        N_Pw = size(pinwheel_tables{base}.x,1);
        pw_ind = zeros([Size_Pw 1]);
        pw_ind(1:N_Pw) = reshape(1: N_Pw,[N_Pw 1]);
    elseif sampleNum ==2
        pw_ind = zeros([Size_Pw 1]);
        ini_ =tracking_tables{2}.ini(:,2);
        end_ =tracking_tables{2}.end(:,2); 
        pw_ind(1:length(ini_)) = ini_;
        pw_ind(length(ini_)+1:length(end_)+length(ini_))=end_;
    else
        LargestLabel = nansum(nansum(pinwheel_stats.label,2)>0,1);
        pw_ind = getIndex(pinwheel_stats.label(:,sampleNum-1),tracking_tables{sampleNum},LargestLabel,Size_Pw);
    end
    
    % fill the table
    pinwheel_stats.label(:,sampleNum) = pw_ind;
    pinwheel_stats.x(pw_ind>0,sampleNum) = pinwheel_tables{sampleNum}.x(pw_ind(pw_ind>0));
    pinwheel_stats.y(pw_ind>0,sampleNum) = pinwheel_tables{sampleNum}.y(pw_ind(pw_ind>0));
    pinwheel_stats.sign(pw_ind>0,sampleNum) = pinwheel_tables{sampleNum}.sign(pw_ind(pw_ind>0));  
end
% pinwheel_stats.probability = sum(pinwheel_stats.label(:,samples)>0,2)/length(samples);%boot_samples;
%   
% % make extra output of pinwheels that are not matched from the other tables
% if nargout == 2
%     pinwheel_spurious = struct('x',[],'y',[],'sign',[],'number',[],'label',[]);
%     
%     for ii=2:length(tracking_tables)
%         if isempty(tracking_tables{ii})
%             continue
%         end
%         
%         ID_table = tracking_tables{ii}(:,ii);
%         ID_pinwheel = pinwheel_tables{ii}.number;
%         
%         [~,test] = intersect(ID_pinwheel,ID_table);
%         
%         pinwheel_spurious.x = [pinwheel_spurious.x;pinwheel_tables{ii}.x(test)];
%         pinwheel_spurious.y = [pinwheel_spurious.y;pinwheel_tables{ii}.y(test)];
%         pinwheel_spurious.sign = [pinwheel_spurious.sign;pinwheel_tables{ii}.sign(test)];
%         pinwheel_spurious.number = [pinwheel_spurious.number;pinwheel_tables{ii}.number(test)];
%         tmp = pinwheel_tables{ii}.number(test);
%         pinwheel_spurious.label = [pinwheel_spurious.label;[ii*ones(size(tmp)) tmp]];
%         
%     end
% end


end

function Ind = getIndex(label,tracking_table,LargestLabel,Size_Pw)
    Ind = zeros([Size_Pw 1]);
    for ii = 1:size(tracking_table.ini,1)
        jj = find(label==tracking_table.ini(ii,1));
        if ~isempty(jj)
            Ind(jj) = tracking_table.ini(ii,2);
        end
    end
    MissingPws =tracking_table.end(:,2); 
    Ind(LargestLabel+1:LargestLabel+length(MissingPws)) = MissingPws;
end


function plotPinwheels(z3D,ROI2D,pinwheel_table,IDs)
    figure;
    plot_map(z3D,ROI2D,0,1)
    hold on
    plot(pinwheel_table.x,pinwheel_table.y,'xw')
    
    label_offset = 1;
    for ii = 1:size(IDs,1)
        ID = IDs(ii);
        if ID ~=0
            x = pinwheel_table.x(ID)+label_offset;
            y = pinwheel_table.y(ID)+label_offset;

            text(x,y,num2str(ii),'Color','white')
        end
    end
end
