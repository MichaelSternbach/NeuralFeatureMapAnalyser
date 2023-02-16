function filter_stats = get_filter_stats(data_obj,tracker_obj,filter_list)
%{
 filter_stats = get_filter_stats(data_obj,tracker_obj,filter_list)

 This code tracks pinwheels from a data source over a range of lowpass
 filter settings. It calculates corresponding pinwheels between the
 filters and extracts their location and sign.

 INPUT
 data_obj: handle to class that does the data handling, eg. extracting
       bootstrap samples, reading ROI, filtering, etc.
       => see class: data_handle
 tracker_obj: handle to class that does the tracking between input maps
       => see class: pinwheel_tracker
 filter_list: range over which pinwheels are tracked. The setting defined
        in the data_obj is included in this range in this code.
        def = linspace(0.2,0.8,50);

 OUTPUT:
 filter_stats: structure with results containing
    - llp_cutoffs : lowpass settings used
    - llp_idx : index to the setting defined in data_obj
    - tracking_table : complex table with corresponding pinwheels
    - pinwheels : x-y-sign information for pinwheels
  
 29/08/2017 : chepe@nld.ds.mpg.de
%}
%%
% read input
if nargin<3
    filter_list = linspace(0.2,0.8,50);
end
llp_cutoffs = unique([data_obj.filter_parameters.lowpass,filter_list]);
llp_idx = find(llp_cutoffs==data_obj.filter_parameters.lowpass);

% get base map
map_raw = data_obj.read_map(1);

% reset tracking
tracker_obj.restart_tracker
tracker_obj.set_parameter('track_annihilated',true)
tracker_obj.set_parameter('track_lost',true)
tracker_obj.set_parameter('add_interpolated_step',true)

% do tracking across cutoffs
pinwheels = struct('x',[],'y',[],'sign',[]);
lp_original = data_obj.filter_parameters.lowpass;
for llp = 1:length(llp_cutoffs)
    % set lowpass
    data_obj.set_filter_parameters('lowpass',llp_cutoffs(llp));
    % track pinwheels
    pw = tracker_obj.add_step(data_obj.filter_map(map_raw),data_obj.ROI);
    % add to lists
    pinwheels(llp).x = pw.x;
    pinwheels(llp).y = pw.y;
    pinwheels(llp).sign = pw.sign;
end
data_obj.set_filter_parameters('lowpass',lp_original);

% prepare output
filter_stats = [];
filter_stats.llp_cutoffs = llp_cutoffs;
filter_stats.llp_idx = llp_idx;
filter_stats.tracking_table = tracker_obj.tracking_table;
filter_stats.pinwheels = pinwheels;

end