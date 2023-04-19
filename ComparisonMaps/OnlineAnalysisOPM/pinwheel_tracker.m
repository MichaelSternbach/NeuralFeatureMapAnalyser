classdef pinwheel_tracker < handle
    % Class to track pinwheels. When initiated a tracking table is started
    % with a given pinwheel list. Then new pinwheel lists can be given as
    % input and the table is expanded. At the end the table can be
    % compressed to match initial and final step. Only requires the x and y
    % coordinates of the pinwheels (if function find_pinwheels is used,
    % those are already given in a structure + assigned pinwheel number).
    %
    % HOW TO USE:
    % Initiate:
    %           obj = pinwheel_tracker()
    % Add map to tracking table:
    %           pw_list = obj.add_step(z,ROI)
    %                   z = complex orientation map, ROI = region of interest
    % Interpolate between two maps in N steps
    %           obj.interpolate(z_ini,z_end,ROI)
    % See full tracking table:
    %           obj.tracking_table
    % Get initial_to_final compressed table
    %           tracking = obj.reduce_tracking
    %
    % OUTPUT:
    %   -obj.tracking_table
    %       column: step
    %       real positive values = pinwheel label in step
    %       real negative values = corresponding label in generated pairs
    %       imaginary negative values = corresponding label in annihilated pairs
    %       zero = lost pinwheel
    %   -tracking:
    %       structure with fields [ini, end]. Each element is a
    %       two column vector showing matching pinwheel numbers. In 'ini' all the
    %       matching, annihilated and lost pinwheels are given. In 'end' all
    %       the pinwheels not related to 'ini' are listed.
    %       Negative signs show the corresponding annihilated/generated pair.
    %
    % PARAMETERS: set with obj.set_parameter(parameter_name,parameter_value,paired_value)
    %  -Max distances for pinwheel movement in a step and when to remove close pinwheels:
    %    limit_movement = 20
    %    limit_interaction = 20
    %    limit_distance = 3
    %  -What to do with lost/annihilated pinwheels or sign changes
    %    track_annihilated = false
    %    track_lost = true
    %    control_sign_change = true
    %    pre_track_events = true
    %  -How to calculate distances
    %    periodic_boundary = false
    %    box_size               Size of the containing box.
    %  -Number of steps when interpolating between two maps
    %    interpolation_steps = 40
    %    add_interpolated_step = false
    %  -Save data of the tracking for debugging?
    %    debug_mode = 0         0: None, 1: Display possible errors
    %                           2: Save all maps for display
    %    debug_dist_ROI           Distance of each pixel to ROI boundary
    %    debug_frames           All maps are saved
    %    debug_limit           When to rise an alarm about missing pinwheels
    %
    % METHODS:
    %   obj.restart_tracker(): empty all tables and start new tracking
    %   pw_list = obj.add_step(z,ROI)       : track a new pinwheel list from given
    %           map and add to the table, returns current pinwheel list
    %   tracking_output = interpolate(z_ini,z_end,ROI) : track pinwheels
    %           while interpolating between the maps. Returns the reduced
    %           tracking table.
    %   pw_list = find_pinwheels(z,ROI,step)  : find pinwheels in given map
    %   [c,a,l,g,n]=obj.match_pinwheels(pw1,pw2): match two pinwheel lists
    %       using the parameters defined in the class (see method description)
    %   tracking = obj.reduce_tracking: reduce table to structure matching
    %       initial and final list
    %   obj.debug_gui: show the tracking steps in GUI to check for errors
    %
    % PRIVATE FUNCTIONS:
    %  -Structure handling for pinwheel lists
    %    structure_create/_delete/_extract/_join/_find
    %  -detect line crossing for pinwheel finder
    %    [crossing,m,n] = line_line_intersection(P1x,P1y,P2x,P2y,Q1x,Q1y,Q2x,Q2y)
    %  -Determine wich compared pinwheels pair belong to each other
    %    adjacency = test_adjacency(pw_distance,possible_adjacency,type)
    %    blocks = separate_adjacency_array(possible_adjacency,type)
    %    model = construct_model(possible_adjacency,costs_array,type)
    %
    % chepe@nld.ds.mg.de 3/Nov/2015
    %%
    
    properties (SetAccess = private)
        % a partial copy of the object
        helper = [];
        
        % Pinwheel lists and base map:
        z_base
        pw_base
        pw_lost
        pw_annihilated
        
        % tracking table and current step
        tracking_table
        step = 0
        
        % reference size for distances
        ref_size
        
        % distance to match moving pinwheels, interacting pairs and when to
        % remove pinwheels that are too close to each other
        limit_movement = 2;%20
        limit_interaction = 2;%20
        
        % pinwheel finder options
        use_external_finder = false
        limit_distance = 2  % use NaN if no pinwheels should be removed
        
        % periodic boundary conditions? Size of the box
        periodic_boundary = false
        box_size
        
        % what kind of pinwheels are tracked after they are not found in
        % one of the steps? Are faulty sign changes taken into
        % consideration? Are possible annihilation/generation removed
        % before the matching?
        pre_track_events = true;
        track_annihilated = false % annih.and gen. are only matched in pairs
        track_lost = true
        control_sign_change = true
        
        % how many steps are used to interpolate between two maps.
        % Are two maps interpolated when adding a step?
        interpolation_steps = 60;%40
        add_interpolated_step = false
        
        % debug mode to plot the tracking process. Needs a distance from
        % ROI measure for every pixel position
        debug_mode = 0
        debug_dist_ROI
        debug_frames = [];
        debug_limit = 40;
    end
    
    methods
        %% ------------- constructor method
        function obj = pinwheel_tracker(ref_size)
            % Constructs a class with pre-defined parameters
            % If started with input argument, the parameters are adjusted
            % to it. ref_size is the typical size of an hypercolumn in
            % pixels.
            
            % initiate step counting
            obj.step = 0;
            
            % pinwhell list layout
            layout = struct('x',[],'y',[],'sign',[],'number',[],'label',[]);
            
            % create new structures for future pinwheels that are not matched
            obj.pw_lost = structure_create(layout);
            obj.pw_annihilated = structure_create(layout,'partner');
            
            % match distances to ref_size
            if nargin==1
                obj.ref_size = ref_size;
                obj.limit_movement = 0.25*obj.ref_size;
                obj.limit_interaction = 0.25*obj.ref_size;
                obj.debug_limit = 0.5*obj.ref_size;
            end
        end
        
        function start_helper(obj)
            
            % delete previous helper
            obj.delete_helper
            
            % Instantiate new object of the same class.
            obj.helper = feval(class(obj));
            
            % Copy all non-hidden properties.
            p = properties(obj);
            for i = 1:length(p)
                obj.helper.(p{i}) = obj.(p{i});
            end
            
            disp('A copy of the class was generated')
        end
        
        function delete_helper(obj)
            obj.helper = [];
        end
        
        function set_parameter(obj,parameter_name,parameter_value,varargin)
            if obj.step>1
                warning('Tracking already started!')
            end
            switch parameter_name
                case 'ref_size'
                    obj.ref_size = parameter_value;
                    obj.limit_movement = 0.25*obj.ref_size;
                    obj.limit_interaction = 0.25*obj.ref_size;
                    obj.debug_limit = 0.5*obj.ref_size;
                case 'limit_movement'
                    obj.limit_movement = parameter_value;
                case 'limit_interaction'
                    obj.limit_interaction = parameter_value;
                case 'use_external_finder'
                    obj.use_external_finder = parameter_value;
                case 'limit_distance'
                    obj.limit_distance = parameter_value;
                case 'periodic_boundary'
                    if isempty(varargin)
                        error('Give the ROI as second input')
                    end
                    obj.periodic_boundary = parameter_value;
                    if all(size(varargin{1})==[1 2])
                        obj.box_size = varargin{1};
                    else
                        obj.box_size = size(varargin{1});
                    end
                case 'track_annihilated'
                    obj.track_annihilated = parameter_value;
                case 'track_lost'
                    obj.track_lost = parameter_value;
                case 'pre_track_events'
                    obj.pre_track_events = parameter_value;
                case 'add_interpolated_step'
                    obj.add_interpolated_step = parameter_value;
                case 'interpolation_steps'
                    obj.interpolation_steps = parameter_value;
                case 'debug_mode'
                    if isempty(varargin)
                        error('Give the ROI as second input')
                    end
                    obj.debug_mode = parameter_value;
                    obj.debug_dist_ROI = bwdist(~varargin{1});
                case 'debug_limit'
                    obj.debug_limit = parameter_value;
                otherwise
                    disp('Unrecognized paramater to set!')
            end
        end
        
        function restart_tracker(obj,complete_restart)
            % restart tracker with already defined parameters
            
            if nargin==1
                complete_restart = true;
            end
            
            % create new structures for future pinwheels that are not matched
            layout = struct('x',[],'y',[],'sign',[],'number',[],'label',[]);
            obj.pw_lost = structure_create(layout);
            obj.pw_annihilated = structure_create(layout,'partner');
            
            if complete_restart
                
                % clear tables
                obj.step = 0;
                obj.tracking_table = [];
                
                obj.pw_base = [];
                obj.debug_frames = [];
            end
        end
        
        %% ------------- feed new list of pinwheels
        function pw_input = add_step(obj,z_input,ROI,do_interpolation)
            % Code that takes two pinwheel lists samples (pw_base/pw_input) and finds the
            % corresponding pairs. In the comparisson an array of
            % matching pinwheels between the interpolations is calculated. A structure of
            % missing pinwheels during the steps (pw_lost) is updated to compare
            % it with the new pinwheels found in further steps
            % INPUT:
            % z = complex orientation map
            % ROI = region of interest
            % do_interpolation = are the layouts interpolated when a step
            % is done?
            % OUTPUT:
            % pw_input = pinwheel list of the current orientation map
            % 15/05/2014 -- chepe@nld.ds.mpg.de
            %%
            
            if ~exist('ROI','var')
                ROI = true(size(z_input));
            else
                if isempty(ROI)
                    ROI = true(size(z_input));
                end
            end
            
            if ~exist('do_interpolation','var')
                do_interpolation = obj.add_interpolated_step;
            end
            
            %% 1) get list of pinwheels and start tracking if first step
            
            obj.step = obj.step + 1;
            pw_input = obj.find_pinwheels(z_input,ROI);
            
            if obj.debug_mode > 0
                obj.debug_frames{obj.step,1}=pw_input;
            end
            if obj.debug_mode == 2
                obj.debug_frames{obj.step,2}=z_input;
                obj.debug_frames{obj.step,2}(~ROI)=NaN;
            end
            
            % If first step, start tables
            if obj.step == 1
                obj.tracking_table = zeros(length(pw_input.number),1);% 'int16';
                obj.tracking_table(:,1) = pw_input.number;
                obj.z_base = z_input;
                obj.pw_base = pw_input;
                return
            else
                obj.tracking_table = [obj.tracking_table,zeros(size(obj.tracking_table,1),1)];
            end
            
            %% 2) match pinwheels with base list
            
            % bypass the interpolation if called internally in add step
            % (input parameter 'do_interpolation' =  false)
            if do_interpolation
                
                % start helper class if needed
                if isempty(obj.helper)
                    obj.start_helper
                end
                
                [~,~,tracking_structure] = obj.helper.interpolate(obj.z_base,z_input,ROI);
                
                % convert output:
                c = tracking_structure.matched;
                a = [tracking_structure.annihilated;tracking_structure.annihilated(:,[2 1])];
                g = [tracking_structure.generated;tracking_structure.generated(:,[2 1])];
                l = tracking_structure.lost;
                n = tracking_structure.new;
                
            else
                [c,a,l,g,n]=obj.match_pinwheels(obj.pw_base,pw_input);
            end
            
            % Check lost and new
            if obj.debug_mode > 0
                % Lost pinwheels
                if ~isempty(l)
                    dist = ceil(obj.debug_dist_ROI(sub2ind( size(obj.debug_dist_ROI), ceil(obj.pw_base.y(l)) ,ceil(obj.pw_base.x(l)) )))';
                    idx2=find(dist>obj.debug_limit);
                    if ~isempty(idx2)
                        disp(['Interpolation step ',num2str(obj.step),': lost pinwheels ',mat2str(l(idx2)),' are ',mat2str(dist(idx2)),' pixels away from ROI boundary.'])
                    end
                end
                % New pinwheels
                if ~isempty(n)
                    dist = ceil(obj.debug_dist_ROI(sub2ind( size(obj.debug_dist_ROI), ceil(pw_input.y(n)) ,ceil(pw_input.x(n)) )))';
                    idx2=find(dist>obj.debug_limit);
                    if ~isempty(idx2)
                        disp(['Interpolation step ',num2str(obj.step),': new pinwheels ',mat2str(n(idx2)),' are ',mat2str(dist(idx2)),' pixels away from ROI boundary.'])
                    end
                end
            end
            
            % fill corresponding pinwheels in tracking array
            [~,ind_tracking,ind_c] = intersect(obj.tracking_table(:,obj.step-1),c(:,1));
            obj.tracking_table(ind_tracking,obj.step) = c(ind_c,2);
            
            % mark annihilated pairs to list
            if ~isempty(a)
                [~,ind_tracking,ind_a] = intersect(obj.tracking_table(:,obj.step-1),a(:,1));
                obj.tracking_table(ind_tracking,obj.step) = complex(real(obj.tracking_table(ind_tracking,obj.step)),-a(ind_a,2));
            end
            
            %% 3) match generated pinwheels with previously lost ones
            % Only aplies when BOTH partners are matched!
            if ~isempty(g)
                g_backUp = g;
                
                if obj.track_annihilated
                    
                    % extract generated list and check matching
                    pw_generated=structure_extract(pw_input,g(:,2));
                    c_annihilated_generated=obj.match_pinwheels(obj.pw_annihilated,pw_generated);
                    
                    % add matching to tracking list and delete from new list
                    if ~isempty(c_annihilated_generated)
                        % sort such that annihilated pinwheels to delete are in the right order
                        c_annihilated_generated = sortrows(c_annihilated_generated,-1);
                        
                        % check which pairs satisfy the 'both partners' condition
                        ann_number = obj.pw_annihilated.number(c_annihilated_generated(:,1));
                        ann_label = obj.pw_annihilated.label(c_annihilated_generated(:,1));
                        ann_partner = obj.pw_annihilated.partner(c_annihilated_generated(:,1));
                        
                        gen_number = pw_generated.number(c_annihilated_generated(:,2));
                        gen_label = obj.step;
                        gen_partner = g(c_annihilated_generated(:,2),1);
                        
                        for ii = 1:size(c_annihilated_generated,1)
                            
                            % check that the partner is also in the matched list
                            co_ann = intersect(ann_number,ann_partner(ii));
                            co_gen = intersect(gen_number,gen_partner(ii));
                            
                            if  isempty(co_ann) || isempty(co_gen)
                                continue
                            end
                            if  isempty(intersect([co_ann co_gen],[ann_number gen_number],'rows'))
                                continue
                            end
                            
                            % add matching pinwheel to tracking list
                            ind = find(obj.tracking_table(:,ann_label(ii))==ann_number(ii));
                            obj.tracking_table(ind,gen_label)=gen_number(ii);
                            
                            % remove from new  list and annihilated list
                            g(g(:,2)==gen_number(ii),:)=[];
                            obj.pw_annihilated = structure_delete(obj.pw_annihilated,c_annihilated_generated(ii,1));
                        end
                    end
                end
                
                % add remaining generated pinwheels to tracking list
                if ~isempty(g)
                    obj.tracking_table(end+(1:size(g,1)),obj.step)=g(:,2);
                end
                
                % mark generated pairs to list
                [~,ind_tracking,ind_g] = intersect(obj.tracking_table(:,obj.step),g_backUp(:,2));
                obj.tracking_table(ind_tracking,obj.step-1) = complex(-g_backUp(ind_g,1),imag(obj.tracking_table(ind_tracking,obj.step-1)));
            end
            
            %% 4) match new pinwheels with previously lost ones
            if ~isempty(n)
                
                % extract new list and check matching
                pw_new=structure_extract(pw_input,n);
                
                if obj.track_lost
                    c_lost_new=obj.match_pinwheels(obj.pw_lost,pw_new);
                    
                    % add matching to tracking list and delete from new list
                    if ~isempty(c_lost_new)
                        % sort such that lost pinwheels to delete are in the right order
                        c_lost_new = sortrows(c_lost_new,-1);
                        
                        lost_number = obj.pw_lost.number(c_lost_new(:,1));
                        lost_label = obj.pw_lost.label(c_lost_new(:,1));
                        new_number = pw_new.number(c_lost_new(:,2));
                        new_label = obj.step;
                        
                        for ii = 1:size(c_lost_new,1)
                            % add matching pinwheel to tracking list
                            ind = find(obj.tracking_table(:,lost_label(ii))==lost_number(ii));
                            obj.tracking_table(ind,new_label)=new_number(ii);
                            % remove from new list
                            n(n==new_number(ii))=[];
                        end
                        % remove from lost list
                        obj.pw_lost = structure_delete(obj.pw_lost,c_lost_new(:,1));
                    end
                end
                
                % add remaining new pinwheels to tracking list
                if ~isempty(n)
                    obj.tracking_table(end+(1:length(n)),obj.step)=n;
                end
            end
            
            %% 5) increase lists
            
            % add annihilated pinwheels to pw_annihilated
            if ~isempty(a)
                tmp = structure_extract(obj.pw_base,a(:,1));
                tmp.partner = a(:,2);
                obj.pw_annihilated = structure_join(obj.pw_annihilated,tmp);
            end
            
            % add lost pinwheels to pw_lost
            if ~isempty(l)
                obj.pw_lost = structure_join(obj.pw_lost,structure_extract(obj.pw_base,l(:)));
            end
            
            
            %% 6) shift base pinwheel structure and map
            obj.pw_base = pw_input;
            obj.z_base = z_input;
        end
        
        %% ------------- interpolate between two maps
        function [tracking_out,tracking_table,tracking_structure] = interpolate(obj,z_ini,z_end,ROI)
            % Take two maps, ini and end, and interpolate between them. In
            % each step track the pinwheels
            
            if ~exist('ROI','var')
                ROI = true(size(z_ini));
            end
            
            % restart the tracker
            obj.restart_tracker;
            
            % interpolate between input maps
            interpolation_weight = linspace(0,1,obj.interpolation_steps);
            for int_step = 1:obj.interpolation_steps
                if int_step == 1
                    z_new = z_ini;
                elseif int_step<obj.interpolation_steps
                    z_new = (1-interpolation_weight(int_step))*z_ini + interpolation_weight(int_step)*z_end;
                else
                    z_new = z_end;
                end
                % track the pinwheels without EXTRA interpolation
                obj.add_step(z_new,ROI,false);
            end
            
            % returned the reduce tracking tables
            [tracking_out,tracking_structure] = obj.reduce_tracking;
            tracking_table = obj.tracking_table;
        end
        
        %% ------------- find pinwheels
        function pinwheels = find_pinwheels(obj,z,ROI,varargin)
            % Code that calculates the pinwheel positions in map
            % Two methods are available: 
            % use_external_finder = true : standard method used in the lab
            %                       false: internal method
            % Input:
            % z: 2D complex map
            % ROI: region of interest array, same dimension as z filled with 0/1 (logical)
            % pw_label: numeric label of the map used to find pinwheels (used in pinwheel tracking)
            %
            % Output:
            % pinwheels: structure with elements [total, position_x, position_y, sign, number, label]
            %
            % 23/04/2014 -- chepe@nld.ds.mpg.de
            if nargin==2
               ROI = true(size(z)); 
            end
            
            if isempty(varargin)
                label = obj.step;
            else
                label = varargin{1};
            end
            
            %% in case of periodic boundary conditions, add extra pieces to 
            %% NOT miss pinwheels as they cross box
            
            if obj.periodic_boundary
                % original size
                [size_y,size_x] = size(z);
                sy = floor(size_y/4);
                sx = floor(size_x/4);
                
                % rescaled size
                z = [z(end-sy:end,end-sx:end) ,  z(end-sy:end,:)  ,  z(end-sy:end,1:sx) ;...
                    z(:,end-sx:end)           ,  z                ,  z(:,1:sx) ;...
                    z(1:sy,end-sx:end)        ,  z(1:sy,:)        ,  z(1:sy,1:sx) ];
                
                ROI = [ROI(end-sy:end,end-sx:end) ,  ROI(end-sy:end,:)  ,  ROI(end-sy:end,1:sx) ;...
                      ROI(:,end-sx:end)           ,  ROI                ,  ROI(:,1:sx) ;...
                      ROI(1:sy,end-sx:end)        ,  ROI(1:sy,:)        ,  ROI(1:sy,1:sx) ];
            end
                
            %% External or internal method:
            
            if obj.use_external_finder
                % Use original code
                [~,~,~,PWxList,PWyList,signList] = find_pinwheels(z,0,ROI,false);
                
                %% prepare output
                if isempty(varargin)
                    label = obj.step;
                else
                    label = varargin{1};
                end
                pinwheel_list = [PWxList(:),PWyList(:),signList(:)];
                pinwheel_list=sortrows(pinwheel_list(:,:),1);
                pinwheels.x = pinwheel_list(:,1);
                pinwheels.y = pinwheel_list(:,2);
                pinwheels.sign = sign(pinwheel_list(:,3));
                pinwheels.number = (1:size(pinwheel_list,1))';
                pinwheels.label = label*ones(size(pinwheel_list,1),1);
                
            else
                
                % Step 1: ROI manipulation
                % Step 2: Calculate contour and find maximum and minimum in x and y for each
                % Step 3: Find which contours boundaries are crossing
                % Step 4: Check crossing points to find pinwheels
                
                % remove the warning for the function line_line_intersection
                warning off MATLAB:divideByZero;
                                
                %% STEP 1: convert ROI to NaN
                
                % fill with NaN such that no contour is detected in borders with ROI
                ROI = double(ROI);
                ROI(~ROI)=NaN;
                
                %% STEP 2: calculate contours and organize into arrays
                
                c_imag=contourc(imag(z).*ROI,[0,0]);c_imag_idx = find(c_imag(1,:) == 0);
                c_real=contourc(real(z).*ROI,[0,0]);c_real_idx = find(c_real(1,:) == 0);
                
                % imaginary part
                % format: (min_x, max_x, min_y, max_y) 1xNx4
                bounds_imag = zeros(1,length(c_imag_idx),4);
                for ii_imag = 1:length(c_imag_idx)
                    tmp = c_imag(:,c_imag_idx(ii_imag)+(1:c_imag(2,c_imag_idx(ii_imag))));
                    tmp = [min(tmp,[],2);max(tmp,[],2)];
                    bounds_imag(1,ii_imag,:)  =[tmp(1) tmp(3) tmp(2) tmp(4)];
                end
                
                % real part
                % here a normalization factor is calculated such that when applying to bounds_imag
                % the elements in bounds_real are [min_x,min_y]=-1 and [max_x,max_y]=1
                % format: ( x_min+(xmax-xmin)/2 y_min+(ymax-ymin)/2 (max_x-min_x)/2, (max_y-min_y)/2)  Mx1x4
                bounds_real = zeros(length(c_real_idx),1,4);
                for ii_real = 1:length(c_real_idx)
                    tmp = c_real(:,c_real_idx(ii_real)+(1:c_real(2,c_real_idx(ii_real))));
                    tmp = [min(tmp,[],2);max(tmp,[],2)];
                    bounds_real(ii_real,1,:)  =[tmp(3)+tmp(1) tmp(4)+tmp(2) tmp(3)-tmp(1) tmp(4)-tmp(2)]/2;
                end
                
                %% STEP 3: check which contours are in the same region
                % shift and normalize all elements in bounds_imag by each factor in bounds_real
                % Also, label coordinates that are x<=-1 as -1, -1<x<1 as 0 and x>=1 as 1
                %   MxNx2  =   (  1xNx2   -    Mx1x1  )  /     Mx1x1
                test_x = sign(fix( bsxfun(@rdivide,  bsxfun(@minus,bounds_imag(:,:,[1 2]),bounds_real(:,:,1))  ,bounds_real(:,:,3)) ));
                test_y = sign(fix( bsxfun(@rdivide,  bsxfun(@minus,bounds_imag(:,:,[3 4]),bounds_real(:,:,2))  ,bounds_real(:,:,4)) ));
                
                % intersection when min*max=0 (one of the points are between -1 and 1)
                % or min+max=0 (min<1 and max>1), for x and y.
                c_overlap=(test_x(:,:,1).*test_x(:,:,2)==0 | test_x(:,:,1)+test_x(:,:,2)==0) & (test_y(:,:,1).*test_y(:,:,2)==0 | test_y(:,:,1)+test_y(:,:,2)==0);
                
                %% STEP 4: test the overlapping contours
                
                % prepare pinwheel list variable
                pinwheel_list= [];
                
                % calculate gradient for pinwheel sign
                [grx,gry]=gradient(real(z));
                [gix,giy]=gradient(imag(z));
                
                % [kx,ky]=meshgrid([0:floor(size(z,1)/2)-1,-ceil(size(z,1)/2):-1]);
                % grx = real(ifft2(1i.*kx.*fft2(real(z))));
                % gry = real(ifft2(1i.*ky.*fft2(real(z))));
                % gix = real(ifft2(1i.*kx.*fft2(imag(z))));
                % giy = real(ifft2(1i.*ky.*fft2(imag(z))));
                
                % when two pinwheels are close to each other, rounding can cause the sign
                % of a pinwheel to switch to the other pinwheel. For this, the determinant
                % is ceil'd and floor'd and the highest absolute value is chosen. Then the
                % exact pinwheel position is floor'd
                
                determinant=grx.*giy-gry.*gix;
                % tmp=determinant;
                % tmp_test = circshift(tmp,[-1 0]);
                % determinant(abs(determinant)<abs(tmp_test))=tmp_test(abs(determinant)<abs(tmp_test));
                % tmp_test = circshift(tmp,[0 -1]);
                % determinant(abs(determinant)<abs(tmp_test))=tmp_test(abs(determinant)<abs(tmp_test));
                % tmp_test = circshift(tmp,[-1 -1]);
                % determinant(abs(determinant)<abs(tmp_test))=tmp_test(abs(determinant)<abs(tmp_test));
                
                for test_contour = find(c_overlap(:))'
                    
                    [cont_real,cont_imag] = ind2sub(size(c_overlap),test_contour);
                    
                    coord_real = c_real(:,c_real_idx(cont_real)+(1:c_real(2,c_real_idx(cont_real))))';
                    coord_imag = c_imag(:,c_imag_idx(cont_imag)+(1:c_imag(2,c_imag_idx(cont_imag))))';
                    
                    dist_limit = 2*sqrt(max([diff(coord_real(:,1)).^2 + diff(coord_real(:,2)).^2 ;...
                        diff(coord_imag(:,1)).^2 + diff(coord_imag(:,2)).^2]));
                    
                    %toCheck = find(euclidian_distance(coord_real(1:end-1,:),coord_imag(1:end-1,:))<=dist_limit);
                    toCheck = find(pdist2(coord_real(1:end-1,:),coord_imag(1:end-1,:))<=dist_limit);
                    
                    if ~isempty(toCheck)
                        [list_real,list_imag] = ind2sub([length(coord_real(1:end-1,1)) length(coord_imag(1:end-1,1))],toCheck);
                        
                        [crossing,m] = line_line_intersection(...
                            coord_real(list_real,1),coord_real(list_real,2),...
                            coord_real(list_real+1,1),coord_real(list_real+1,2),...
                            coord_imag(list_imag,1),coord_imag(list_imag,2),...
                            coord_imag(list_imag+1,1),coord_imag(list_imag+1,2));
                        
                        x = coord_imag(list_imag,1)+m.*(coord_imag(list_imag+1,1)-coord_imag(list_imag,1));
                        y = coord_imag(list_imag,2)+m.*(coord_imag(list_imag+1,2)-coord_imag(list_imag,2));
                        
                        pinwheel_list = [pinwheel_list; x(crossing) y(crossing) (determinant(sub2ind(size(determinant),round(y(crossing)),round(x(crossing)))))];
                    end
                    
                end
                
                % remove pinwheels that are in the same pixel (double counts)
                [ix,~]=find( (pdist2(floor(pinwheel_list(:,1)),floor(pinwheel_list(:,1)))+2*eye(size(pinwheel_list,1)))==0 & (pdist2(floor(pinwheel_list(:,2)),floor(pinwheel_list(:,2)))+2*eye(size(pinwheel_list,1)))==0);
                pinwheel_list(ix(1:2:end),:)=[];
                
                %% check if pinwheels are too close to each other
                % (problems when determining their charge)
                
                if ~isnan(obj.limit_distance)
                    
                    % calculate distance
                    % note: periodic conditions already implemented in
                    %       copies of the box on each side
                    X1 = pinwheel_list(:,1);
                    Y1 = pinwheel_list(:,2);                    
                    pw_distance=pdist2([X1 Y1],[X1 Y1],'euclidean');
                    
                    % check pinwheels that are too close
                    possible_adjacency = (pw_distance<obj.limit_distance & eye(length(pinwheel_list(:,1)))~=1);
                    
                    if sum(possible_adjacency(:))>0
                        % eliminate only pairs of pinwheels (the closes ones)
                        adjacency = test_adjacency(pw_distance,possible_adjacency,'symmetric','pw_finder');
                        [~,near_idx] = find(adjacency);
                        pinwheel_list(near_idx,:)=[];
                    end
                end
                
                %% prepare output
                pinwheel_list=sortrows(pinwheel_list(:,:),1);
                pinwheels.x = pinwheel_list(:,1);
                pinwheels.y = pinwheel_list(:,2);
                pinwheels.sign = sign(pinwheel_list(:,3));
                pinwheels.number = (1:size(pinwheel_list,1))';
                pinwheels.label = label*ones(size(pinwheel_list,1),1);
                
            end
                        
            %% in case of periodic boundary conditions, remove extra pieces
            if obj.periodic_boundary
                % original size = size_y, size_x
                % to remove = sy, sx
                
                idx = pinwheels.x < sx | pinwheels.x >= size_x+sx | pinwheels.y < sy | pinwheels.y >= size_y+sy;
                                
                pinwheels.x(idx) = [];
                pinwheels.x = pinwheels.x - sx;
                
                pinwheels.y(idx) = [];
                pinwheels.y = pinwheels.y - sy;
                
                pinwheels.sign(idx) = [];
                
                pinwheels.number = (1:length(pinwheels.x))';
                pinwheels.label(idx) = [];                
            end
            
        end
        
        %% ------------- match pinwheels
        function [c,a,l,g,n]=match_pinwheels(obj,pw1,pw2)
            % Function that gets two lists of pinwheels and finds correspondance
            % depending on the euclidian distance between them and their sign. It also
            % gives a list of annihilated, generated, lost and new pinwheels
            % STEPS:
            % 1) Find possible annihilated/generated pairs (limit interaction/2)
            % 2) Find pre-post corresponding pairs
            % 3) Find annihilated pinwheel pairs
            % 4) Find generated pinwheel pairs
            % Before defining the rest of the pinwheels as lost:
            % 5) Control sign change: corresponding pinwheels (limit movement/2)
            % 6) Control sign change: annihilated pinwheel pairs (limit interaction/2)
            % 7) Control sign change: generated pinwheel pairs (limit interaction/2)
            % And only then
            % 8) Define rest as lost/new pinwheels
            %
            % INPUT:
            % pw1/pw2 = structures of pinwheel data. Used is pw.x, pw.y, pw.sign
            %
            % OUTPUT:
            % c = list of corresponding pinwheels. First column for pw1, second for pw2
            % a = list of annihilated pinwheels pairs in pw1. List in first column gets
            %      annihilated by corresponding list in second column
            % l = list of pinwheels in pw1 that are neither in 'c' nor in 'a'
            % g = list of generated pinwheel pairs in pw2. List in first column gets
            %      generated with corresponding list in second column
            % n = list of pinwheels in pw2 that are neither in 'c' nor in 'g'
            %
            % 25/04/2014 -chepe@nld.ds.mpg.de
            
            % Prepare output variables
            c=[];%int16.empty;
            a=[];%int16.empty;
            l=[];%int16.empty;
            g=[];%int16.empty;
            n=[];%int16.empty;
            
            % Check input parameters
            
            % if any of the pinwheel lists is empty, return empty arrays
            if (isempty(pw1.x)||isempty(pw2.x))
                return
            end
            
            % calculate distances between pinwheels
            % use periodic boundaries if required
            %pw_distance = euclidian_distance([pw1.x pw1.y],[pw1.x pw1.y]);
            X1 = pw1.x;
            Y1 = pw1.y;
            X2 = pw2.x;
            Y2 = pw2.y;
            
            if ~obj.periodic_boundary
                pw_distance_pre=pdist2([X1 Y1],[X1 Y1],'euclidean');
                pw_distance_pos=pdist2([X2 Y2],[X2 Y2],'euclidean');
                pw_distance_crossed=pdist2([X1 Y1],[X2 Y2],'euclidean');
            else                
                dx = abs(bsxfun(@minus,X1,X1'));
                dy = abs(bsxfun(@minus,Y1,Y1'));                
                dx(dx > obj.box_size(2)/2) = obj.box_size(2) - dx(dx > obj.box_size(2)/2);
                dy(dy > obj.box_size(1)/2) = obj.box_size(1) - dy(dy > obj.box_size(1)/2);
                pw_distance_pre = sqrt(dx.^2 + dy.^2);
                
                dx = abs(bsxfun(@minus,X2,X2'));
                dy = abs(bsxfun(@minus,Y2,Y2'));                
                dx(dx > obj.box_size(2)/2) = obj.box_size(2) - dx(dx > obj.box_size(2)/2);
                dy(dy > obj.box_size(1)/2) = obj.box_size(1) - dy(dy > obj.box_size(1)/2);
                pw_distance_pos = sqrt(dx.^2 + dy.^2);
                
                dx = abs(bsxfun(@minus,X1,X2'));
                dy = abs(bsxfun(@minus,Y1,Y2'));                
                dx(dx > obj.box_size(2)/2) = obj.box_size(2) - dx(dx > obj.box_size(2)/2);
                dy(dy > obj.box_size(1)/2) = obj.box_size(1) - dy(dy > obj.box_size(1)/2);
                pw_distance_crossed = sqrt(dx.^2 + dy.^2);               
            end
            
            %% 1) Find possible annihilated/generated pairs (limit interaction/4)
            % --> this is done to get rid of the sign changing problem when
            % the pinwheels get too close. Here the pair is removed in
            % advanced if there is no possible partner after/before the event
            
            if obj.pre_track_events
                
                %----------- find possible annihilated
                possible_adjacency = pw_distance_pre<obj.limit_interaction/4 & bsxfun(@times,pw1.sign,pw1.sign')==-1;
                adjacency = test_adjacency(pw_distance_pre,possible_adjacency,'symmetric','pre_annihilated');
                [ann_1,ann_2]=ind2sub(size(adjacency),find(adjacency));
                
                %----------- find possible generated
                possible_adjacency = pw_distance_pos<obj.limit_interaction/4 & bsxfun(@times,pw2.sign,pw2.sign')==-1;
                adjacency = test_adjacency(pw_distance_pos,possible_adjacency,'symmetric','pre_generated');
                [gen_1,gen_2]=ind2sub(size(adjacency),find(adjacency));
                
                % ------- check when there is no other pinwheel found in the region
                possible_adjacency = pw_distance_crossed<obj.limit_interaction & bsxfun(@times,pw1.sign,pw2.sign')==1;
                
                % annihilated
                ind = sum(possible_adjacency(ann_1,:),2)==0 & sum(possible_adjacency(ann_2,:),2)==0;
                a = [ann_1(ind),ann_2(ind)];%int16([ann_1(ind),ann_2(ind)]);
                
                % generated
                ind = sum(possible_adjacency(:,gen_1),1)==0 & sum(possible_adjacency(:,gen_2),1)==0;
                g=[gen_1(ind),gen_2(ind)];%int16([gen_1(ind),gen_2(ind)]);
                
                % -- Make a list of missing pinwheels
                missing_pw1 = setdiff(1:length(pw1.x),a(:))';
                missing_pw2 = setdiff(1:length(pw2.x),g(:))';
                
            else
                missing_pw1 = [1:length(pw1.x)]';
                missing_pw2 = [1:length(pw2.x)]';
            end
            
            %% 2) Find pre-post corresponding pairs
            
            if ~isempty(missing_pw1) && ~isempty(missing_pw2)
                
                % test which pinwheels with the SAME sign are in the neighborhood and
                % calculate adjacency array
                possible_adjacency = pw_distance_crossed(missing_pw1,missing_pw2)<obj.limit_movement & bsxfun(@times,pw1.sign(missing_pw1),pw2.sign(missing_pw2)')==1;
                adjacency = test_adjacency(pw_distance_crossed(missing_pw1,missing_pw2),possible_adjacency,'asymmetric','matched');
                
                % return indices of found pinwheel pairs
                [c1,c2]=ind2sub(size(adjacency),find(adjacency));
                c = [missing_pw1(c1) , missing_pw2(c2)];%int16([missing_pw1(c1)',missing_pw2(c2)']);
                
                % -- Shorten list of missing pinwheels
                if ~isempty(c)
                    missing_pw1 = setdiff(missing_pw1,c(:,1));
                    missing_pw2 = setdiff(missing_pw2,c(:,2));
                end
                
            end
            
            %% If no more output is required, stop the matching process
            
            if nargout==1
                return;
            end
            
            %% 3) Find annihilated pinwheel pairs
            
            if ~isempty(missing_pw1)
                
                % test which pinwheels with OPPOSITE sign are in the neighborhood and
                % calculate adjacency matrix
                possible_adjacency = pw_distance_pre(missing_pw1,missing_pw1)<obj.limit_interaction & bsxfun(@times,pw1.sign(missing_pw1),pw1.sign(missing_pw1)')==-1;
                adjacency = test_adjacency(pw_distance_pre(missing_pw1,missing_pw1),possible_adjacency,'symmetric','annihilated');
                
                % return found pinwheel indices
                [a1,a2]=ind2sub(size(adjacency),find(adjacency));
                a_tmp = [missing_pw1(a1) , missing_pw1(a2)];
                
                % add to list and remove annihilated pinwheels from missing_pw1 list
                if ~isempty(a_tmp)
                    a = [a;a_tmp];%int16([a;a_tmp]);
                    missing_pw1 = setxor(a_tmp(:,1),missing_pw1);%int16(setxor(a_tmp(:,1),missing_pw1));
                end
                
            end
            
            %% 4) Find generated pinwheel pairs
            
            if ~isempty(missing_pw2)
                
                % test which pw with OPPOSITE sign are in the neighborhood and
                % calculate adjacency matrix
                possible_adjacency = pw_distance_pos(missing_pw2,missing_pw2)<obj.limit_interaction & bsxfun(@times,pw2.sign(missing_pw2),pw2.sign(missing_pw2)')==-1;
                adjacency = test_adjacency(pw_distance_pos(missing_pw2,missing_pw2),possible_adjacency,'symmetric','generated');
                
                % return found pinwheel indices
                [g1,g2]=ind2sub(size(adjacency),find(adjacency));
                g_tmp = [missing_pw2(g1) , missing_pw2(g2)];
                
                % add to list and remove generated pinwheels from missing_pw2 list
                if ~isempty(g_tmp)
                    g = [g;g_tmp];%int16([g;g_tmp]);
                    missing_pw2 = setxor(g_tmp(:,1),missing_pw2);% int16(setxor(g_tmp(:,1),missing_pw2));
                end
                
            end
            
            %%
            
            if obj.control_sign_change
                
                %% 5) Control sign change: corresponding pinwheels (limit movement/2)
                
                if ~isempty(missing_pw1) && ~isempty(missing_pw2)
                    
                    possible_adjacency = pw_distance_crossed(missing_pw1,missing_pw2)<obj.limit_movement/2 & bsxfun(@times,pw1.sign(missing_pw1),pw2.sign(missing_pw2)')==-1;
                    adjacency = test_adjacency(pw_distance_crossed(missing_pw1,missing_pw2),possible_adjacency,'asymmetric','sign_matched');
                    
                    % return indices of found pinwheel pairs
                    [c1,c2]=ind2sub(size(adjacency),find(adjacency));
                    c_tmp = [missing_pw1(c1) , missing_pw2(c2)];
                    
                    % Add to list and shorten list of missing pinwheels
                    if ~isempty(c_tmp)
                        try
                            c = [c;c_tmp];%int16([c;c_tmp]);
                        catch
                            perro = 1;
                        end
                        
                        missing_pw1 = setdiff(missing_pw1,c_tmp(:,1));
                        missing_pw2 = setdiff(missing_pw2,c_tmp(:,2));
                    end
                end
                
                %% 6) Control sign change: annihilated pinwheel pairs (limit interaction/2)
                
                if ~isempty(missing_pw1)
                    
                    % test which pinwheels with OPPOSITE sign are in the neighborhood and
                    % calculate adjacency matrix
                    possible_adjacency = pw_distance_pre(missing_pw1,missing_pw1)<obj.limit_interaction/2 & bsxfun(@times,pw1.sign(missing_pw1),pw1.sign(missing_pw1)')==1 & eye(length(missing_pw1))==0;
                    adjacency = test_adjacency(pw_distance_pre(missing_pw1,missing_pw1),possible_adjacency,'symmetric','sign_annihilated');
                    
                    % return found pinwheel indices
                    [a1,a2]=ind2sub(size(adjacency),find(adjacency));
                    a_tmp = [missing_pw1(a1) ,missing_pw1(a2)];
                    
                    % add to list and remove annihilated pinwheels from missing_pw1 list
                    if ~isempty(a_tmp)
                        a = [a;a_tmp];%int16([a;a_tmp]);
                        missing_pw1 = setxor(a_tmp(:,1),missing_pw1);%int16(setxor(a_tmp(:,1),missing_pw1));
                    end
                    
                end
                
                %% 7) Control sign change: generated pinwheel pairs (limit interaction/2)
                
                if ~isempty(missing_pw2)
                    
                    % test which pw with OPPOSITE sign are in the neighborhood and
                    % calculate adjacency matrix
                    possible_adjacency = pw_distance_pos(missing_pw2,missing_pw2)<obj.limit_interaction/2 & bsxfun(@times,pw2.sign(missing_pw2),pw2.sign(missing_pw2)')==1 & eye(length(missing_pw2))==0;
                    adjacency = test_adjacency(pw_distance_pos(missing_pw2,missing_pw2),possible_adjacency,'symmetric','sign_generated');
                    
                    % return found pinwheel indices
                    [g1,g2]=ind2sub(size(adjacency),find(adjacency));
                    g_tmp = [missing_pw2(g1) , missing_pw2(g2)];
                    
                    % add to list and remove generated pinwheels from missing_pw2 list
                    if ~isempty(g_tmp)
                        g = [g;g_tmp];%int16([g;g_tmp]);
                        missing_pw2 = setxor(g_tmp(:,1),missing_pw2);%int16(setxor(g_tmp(:,1),missing_pw2));
                    end
                    
                end
            end
            %% 8) Define rest as lost/new pinwheels
            
            l =  missing_pw1;
            n =  missing_pw2;
            
            %%
        end
        
        %% ------------- reduce tracking array
        function [tracking_out,pw_result] = reduce_tracking(obj,tracking,return_global_idx)
            % Reduce tracking and fill list of corresponding pinwheels
            % This code reads a tracking array and makes a reduced version of it showing
            % which pinwheels correspond to which (or gets annihilated/generated lost/new)
            % There are 5 possible cases. The pinwheel can:
            % 1) be matched from beginning to end
            % 2) get lost during the interpolations
            % 3) get annihilated
            % 4) appear during the interpolations
            % 5) get generated
            % For the cases 3) and 5), sometimes a pinwheel gets annihilated/generated
            % by a member of a pair that survives until the end/beginning of the
            % interpolation. In this case that pinwheel is tracked. The
            % generation/annihilation of pairs can also happen multiple times during
            % the interpolations, such that this check has to be done recursively
            %
            % If return_global_idx = true:
            % Code returns the indices (rows) of the pinwheels, not the
            % number. If a pinwheel gets lost during annihilation process,
            % the last tracked index is returned (instead of a zero) as
            % negative complex number.
            
            if nargin<2
                tracking = double(obj.tracking_table);                
            end
            
            if nargin<3
                return_global_idx = false;
            end
            
            %% PART 1: Clean tracking table
            
            % make a fake pinwheel to close the loops and crosses in tracking
            % IT should be high enough!!! Larger than the number of pw
            fake_label = max(real(tracking(:)))+1000;
            
            % original number of rows without helpers
            num_rows = size(tracking,1);
            loop_label = fake_label+1;
            
            % ==== clean boarders of tracking table
            % remove annihilation marks before the first step and
            % generation marks after the last step
            tracking(:,1) = tracking(:,1) - 1i*imag(tracking(:,1));
            ind = real(tracking(:,end))<0;
            tracking(ind,end) = tracking(ind,end) - real(tracking(ind,end));
            
            % ====  remove annihilation/generation loops of variable size
            % [ A A -iB 0 0 -B A A ]  -->  [ A A F F F F A A ]
            % [ B B -iA 0 0 -A B B ]       [ B B F F F F B B ]
            
            for loop_length = 0:size(tracking,2)
                
                % find annih-generation in same row
                if loop_length==0
                    [iy,ix]=find(imag(tracking(:,1:end-loop_length))<0  & real(tracking(:,1+loop_length:end))<0 );
                else
                    [iy,ix]=find(imag(tracking(:,1:end-loop_length))<0 & real(tracking(:,1:end-loop_length))==0 & real(tracking(:,1+loop_length:end))<0 & imag(tracking(:,1+loop_length:end))==0);
                    % must have zeros inbetween
                    for ii=length(ix):-1:1
                        if any(tracking(iy(ii),ix(ii)+1:ix(ii)+loop_length-1)~=0)
                            iy(ii) = []; ix(ii) = [];
                        end
                    end
                end
                % jump if there is nothing to do
                if isempty(iy)
                   continue 
                end
                
                % make corresponding ann/gen
                loop = 1i*imag(tracking(sub2ind(size(tracking),iy,ix))) + real(tracking(sub2ind(size(tracking),iy,ix+loop_length)));
                ixm1 = tracking(sub2ind(size(tracking),iy,ix-1));
                ixp1 = tracking(sub2ind(size(tracking),iy,ix+1+loop_length));
                anti_loop = -ixp1 -1i*ixm1;
                
                % get matching loop pairs (have to be simultaneous -> ix)
                [test, idx_c] = ismember([ix loop],[ix anti_loop],'rows');
                % remove
                tmp = 1:length(test);
                loop_pairs = unique(sort([tmp(test)',idx_c(test)],2),'rows');
                for ii=1:size(loop_pairs,1)
                    tracking(iy(loop_pairs(ii,1)),ix(loop_pairs(ii,1))+[0:loop_length]) = fake_label;
                    tracking(iy(loop_pairs(ii,2)),ix(loop_pairs(ii,2))+[0:loop_length]) = fake_label;
                end
            end
            
            % ==== remove pinwheel crosses
            % Cases like this
            % [ A A -C-iB A A ]  ->  [ A A F A A ]
            % [ B B  -iA  0 0 ]      [ B B 0 C C ]
            % [ 0 0  -A   C C ]      [helper B->C]
            
            % find annih/generation in same step
            [iy,ix]=find(real(tracking)<0 & imag(tracking)<0);
            cross = tracking(sub2ind(size(tracking),iy,ix));
            
            % find corresponding ann/gen
            ixm1 = tracking(sub2ind(size(tracking),iy,ix-1));
            ixp1 = tracking(sub2ind(size(tracking),iy,ix+1));
            anti_cross = -ixp1 -1i*ixm1;
            
            % remove matching crosses (should have been removed in previous block)
            test = ismember(cross,anti_cross);
            tracking(sub2ind(size(tracking),iy(test),ix(test))) = fake_label*ones(sum(test),1);
            iy(test) = [];
            ix(test) = [];
            cross(test) = [];
            
            % if there is a cross that happens close to the limits of the
            % tracking, expand such that the matching part to add fits
            if any(ix==size(tracking,2)-1)
                tmp = tracking(:,end);
                tmp(tmp<1)=0;
                tracking = [tracking,tmp];
            end
            if any(ix==2)
                tmp = tracking(:,1);
                tmp(tmp<1)=0;
                tracking = [tmp,tracking];
                ix = ix+1;
            end
            
            % bind ends by adding lines to table
            for looseEnd = 1:length(cross)
                
                % make helper loop
                toAdd = zeros(3,size(tracking,2));
                toAdd(1,ix(looseEnd)+[-2:0]) = [-loop_label-1, loop_label 1i*imag(cross(looseEnd))];
                toAdd(2,ix(looseEnd)+[-2:2]) = [-loop_label, loop_label+1, loop_label+1, loop_label+1, -1i*loop_label];
                toAdd(3,ix(looseEnd)+[0:2]) = [real(cross(looseEnd)), loop_label, -1i*(loop_label+1)];
                
                % find position of loose ends
                iyan = find(tracking(:,ix(looseEnd)-1)==-imag(cross(looseEnd)));
                iyge = find(tracking(:,ix(looseEnd)+1)==-real(cross(looseEnd)));
                
                % do corrections to tracking and add helper loop
                tracking(iy(looseEnd),ix(looseEnd))=fake_label;
                tracking(iyan,ix(looseEnd)) = -1i* loop_label;
                tracking(iyge,ix(looseEnd)) = -loop_label;
                tracking = [tracking;toAdd];
                
                loop_label = loop_label +2;
            end
            
            % ====  remove lost-new loops
            % remove paths of the type
            % [A A 0 0 A A]  ->  [A A F F A A]
            
            num_cols = size(tracking,2);
            tracking(:,num_cols+1) = NaN;
            tracking = tracking.';
            for spurious_length = 1:num_cols-2
                idx =   strfind(sign(real(tracking(:))+imag(tracking(:))).',[1 zeros(1,spurious_length) 1]);
                for toFill=1:length(idx)
                    tracking(idx(toFill)+[1:spurious_length]) = fake_label;
                end
            end
            tracking = tracking.';
            tracking(:,end) = [];
            
            %% PART 2: Go through tracking table row by row
            
            % Start empty array and go row by row
            tracking_out = [];
            try
                for row=1:num_rows
                    
                    % CASE 0: Pinwheel appearing only during interpolation steps
                    % [ 0 % % % % 0 ]
                    
                    if real(tracking(row,1))<1 && real(tracking(row,end))<1
                        continue
                    end
                    
                    % CASE 1: Pinwheel is matched from beginning to end
                    % [ A A A A A A ]
                    
                    if all(real(tracking(row,:))>0)
                        tracking_out = [tracking_out;row row];
                        continue
                    end
                    
                    % CASE 2: Pinwheel gets lost in the same row
                    % [ A A 0 0 0 0 ] <- must be only zeros
                    
                    if real(tracking(row,1))>0 && real(tracking(row,end))==0 && ...
                            all(tracking(row,find(real(tracking(row,:))<1,1,'first'):end)==0)
                        tracking_out = [tracking_out;row -1i*row];
                        continue
                    end
                    
                    % CASE 3: Pinwheel appears in the same row
                    % [ 0 0 0 0 A A ] <- must be only zeros
                    
                    if real(tracking(row,1))==0 && real(tracking(row,end))>0 && ...
                            all(tracking(row,1:find(real(tracking(row,:))<0,1,'last'))==0)
                        tracking_out = [tracking_out;-1i*row row];
                        continue
                    end
                    
                    % CASE 4: Track pinwheel from the beginning
                    % [ A A -iB % % % ] <- follow loop to B
                    
                    if real(tracking(row,1))>0
                        
                        % initiate annihilation - generation loop
                        row_loop = row;
                        continue_loop = true;
                        loop_depth = 0;
                        test_when = 1;
                        while continue_loop
                            
                            if loop_depth==10
                                disp('WARNING: more than 10 generation/annihilation loops!');
                            end
                            
                            % find co-annihilated
                            test_when = find(imag(tracking(row_loop,test_when:end))<0,1,'first')+test_when-1;
                            test_with = find(imag(tracking(:,test_when))==-real(tracking(row_loop,test_when-1)) & real(tracking(:,test_when-1))>0);
                            
                            % check if the co-annihilated has another event
                            if all(real(tracking(test_with,1:test_when-1))>=0 & imag(tracking(test_with,1:test_when-1))==0)
                                % check if it is lost or an annihilation
                                if real(tracking(test_with,1))>0
                                    tracking_out = [tracking_out;row -test_with];
                                else
                                    tracking_out = [tracking_out;row -1i*test_with];
                                end
                                continue_loop = false;
                                
                            else % there is a generation event before the start
                                
                                % find co-generated
                                test_when = find(real(tracking(test_with,1:test_when-1))<0,1,'last');
                                test_with = find(real(tracking(:,test_when))==-real(tracking(test_with,test_when+1)) & real(tracking(:,test_when+1))>0);
                                
                                % check if the co-generated partner has another event
                                if all(real(tracking(test_with,test_when+1:end))>=0 & imag(tracking(test_with,test_when+1:end))==0)
                                    % check if it is lost or a match
                                    if real(tracking(test_with,end))>0
                                        tracking_out = [tracking_out;row test_with];
                                    else
                                        tracking_out = [tracking_out;row -1i*test_with];
                                    end
                                    continue_loop = false;
                                    
                                else  % continue loop with new pinwheel
                                    row_loop = test_with;
                                    loop_depth = loop_depth + 1;
                                    continue_loop = true;
                                end
                                
                            end
                        end
                    end
                    
                    % CASE 5: Track pinwheel from the end
                    % [ % % -iB A A ] <- follow loop to B
                    
                    if real(tracking(row,end))>0
                        
                        % initiate generation - annihilation loop
                        row_loop = row;
                        continue_loop = true;
                        loop_depth = 0;
                        test_when = size(tracking,2);
                        while continue_loop
                            
                            if loop_depth == 10
                                disp('WARNING: more than 10 generation/annihilation loops!');
                            end
                            
                            % find co-generated
                            test_when = find(real(tracking(row_loop,1:test_when))<0,1,'last');
                            test_with = find(real(tracking(:,test_when))==-real(tracking(row_loop,test_when+1)) & real(tracking(:,test_when+1))>0);
                            
                            % check if the co-generated partner exists until the end
                            if all(real(tracking(test_with,test_when+1:end))>=0 & imag(tracking(test_with,test_when+1:end))==0)
                                % check if it is new or a generation
                                if real(tracking(test_with,end))>0
                                    tracking_out = [tracking_out;-test_with row];
                                else
                                    tracking_out = [tracking_out;-1i*test_with row];
                                end
                                continue_loop = false;
                                
                            else % there is an annihilation event before the end
                                
                                % find co-annihilated
                                test_when = find(imag(tracking(test_with,test_when+1:end))<0,1,'first')+test_when;
                                test_with = find(imag(tracking(:,test_when))==-real(tracking(test_with,test_when-1)) & real(tracking(:,test_when-1))>0);
                                
                                % check if the co-annihilated partner has another event
                                if all(real(tracking(test_with,1:test_when-1))>=0 & imag(tracking(test_with,1:test_when-1))==0)
                                    % check if it is new or a match
                                    if real(tracking(test_with,1))>0
                                        tracking_out = [tracking_out;test_with row];
                                    else
                                        tracking_out = [tracking_out;-1i*test_with row];
                                    end
                                    continue_loop = false;
                                    
                                else % continue loop with new pinwheel
                                    row_loop = test_with;
                                    loop_depth = loop_depth + 1;
                                    continue_loop = true;
                                end
                                
                            end
                        end
                    end
                    
                end % END of loop over CASES in the row
            catch err
                % problem when reducing the tracking array
                disp('Error when simplifying tracking array!')
                rethrow(err)
            end
            
            % Remove double tracking from going back and forth
            tracking_out = unique(tracking_out,'rows');
            
            %% PART 3: Get matching directly by row number and compare with previous
            % In making the tracking table, lost and annihilated pinwheels
            % can be re-assigned to new and generated ones depending on the
            % last position they had. This is important for pinwheels that
            % annihilate-generate multiple times in the data. 
            % ==> Take into account rows (=pw) that are matched in a direct
            % reading, but are NOT following the tracking table. In those
            % cases the position-tracking could be wrong (it shouldn't be
            % matched) or incomplete (the other matching parter failed)
            
            % list of corrections types to apply (see below)
            correct_LN = true;
            correct_2LG = true;
            correct_LG = false;
            correct_AG = true;
            correct_MG = false;
            
            % make separate list of results for easier manipulation
            pw_result = [];
            
            pw_result.matched = tracking_out(real(tracking_out(:,1))>0 & real(tracking_out(:,2))>0,:);
            pw_result.lost = tracking_out(real(tracking_out(:,2))==0,1);
            pw_result.new = tracking_out(real(tracking_out(:,1))==0,2);
            
            tmp = real(tracking_out(:,1))>0 & real(tracking_out(:,2))<0;
            pw_result.annihilated = unique(sort(abs(tracking_out(tmp,:)),2),'rows');
            
            tmp = real(tracking_out(:,1))<0 & real(tracking_out(:,2))>0;
            pw_result.generated = unique(sort(abs(tracking_out(tmp,:)),2),'rows');
            
            % get alternative matching based only on the row number
            alternative_matching = find(all(real(tracking(:,[1 end]))>0,2));
            alternative_matching(ismember([alternative_matching alternative_matching],pw_result.matched,'rows'),:) = [];
            
            % 2LG: Two pinwheels get lost and then generated (happens when
            %           their direct connection was not identified)
            % [ P P 0 -A P ] -?-> [ P P 0 0 P P ]
            % [ A A 0 -P A ]      [ A A 0 0 A A ]
            
            to_match = pw_result.generated(all(ismember(pw_result.generated,pw_result.lost),2),:);
            if ~isempty(to_match)
                if correct_2LG
                    pw_result.lost(ismember(pw_result.lost,to_match(:)))=[];
                    pw_result.generated(all(ismember(pw_result.generated,to_match),2),:)=[];
                    pw_result.matched = [pw_result.matched; to_match(:),to_match(:)];
                end
                alternative_matching(ismember(alternative_matching,to_match)) = [];
            end
            
            to_match = pw_result.annihilated(all(ismember(pw_result.annihilated,pw_result.new),2),:);
            if ~isempty(to_match)
                if correct_2LG
                    pw_result.new(ismember(pw_result.new,to_match(:)))=[];
                    pw_result.annihilated(all(ismember(pw_result.annihilated,to_match),2),:)=[];
                    pw_result.matched = [pw_result.matched; to_match(:),to_match(:)];
                end
                alternative_matching(ismember(alternative_matching,to_match)) = [];
            end
            
            % go over remaining alternative matchings
            for aii=length(alternative_matching):-1:1
                row = alternative_matching(aii);
                
                % LN: Same lost and new pinwheel (happens when there are
                %       annihilation/generation in the tracking but those
                %       get lost)
                % [ P P -B 0 0 -C P P ]       [ P P 0 0 0 0 P P ]
                % [ 0 B -P 0 0  0 0 0 ] -?->  [ 0 0 0 0 0 0 0 0 ]
                % [ 0 0  0 0 0 -P C 0 ]       [ 0 0 0 0 0 0 0 0 ]
                
                if any(pw_result.lost==row) && any(pw_result.new==row)
                    if correct_LN
                        [to_match,idxl,idxn] = intersect(pw_result.lost,pw_result.new);
                        pw_result.matched = [pw_result.matched; to_match,to_match];
                        pw_result.lost(idxl) = [];
                        pw_result.new(idxn) = [];
                    end
                    alternative_matching(aii) = []; continue
                end
                
                % LG: Lost and generated (happens close to ROI border)
                % [ P P 0 -A P] -?-> [ P P 0 0 P ]
                % [ 0 0 0 -P A]      [ 0 0 0 0 A ]
                
                if any(pw_result.lost==row) && any(pw_result.generated(:)==row)
                    if correct_LG
                        pw_result.lost(pw_result.lost == row) = [];
                        
                        ind = any(pw_result.generated==row,2);
                        pw_result.new = [pw_result.new;pw_result.generated(ind,pw_result.generated(ind,:)~=row)];
                        pw_result.generated(ind,:) = [];
                        
                        pw_result.matched = [pw_result.matched; row,row];
                    end
                    alternative_matching(aii) = []; continue
                    
                elseif any(pw_result.new==row) && any(pw_result.annihilated(:)==row)
                    if correct_LG
                        pw_result.new(pw_result.new == row) = [];
                        
                        ind = any(pw_result.annihilated==row,2);
                        pw_result.lost = [pw_result.lost;pw_result.annihilated(ind,pw_result.annihilated(ind,:)~=row)];
                        pw_result.annihilated(ind,:) = [];
                        
                        pw_result.matched = [pw_result.matched; row,row];
                    end
                    alternative_matching(aii) = []; continue
                end
                
                % AG: Annihilated and generated (happens when the new pair
                %                       is too far away from the old one)
                % [ A -P 0  0 0 ]       [ A -P 0 -P B ]
                % [ P -A 0 -B P ] -?->  [ P -A 0 -B P ]
                % [ 0  0 0 -P B ]       [ 0  0 0  0 0 ]
                
                if any(pw_result.annihilated(:)==row) && any(pw_result.generated(:)==row)
                    if correct_AG
                        ind1 = any(pw_result.annihilated==row,2);
                        ind2 = any(pw_result.generated==row,2);
                        
                        pw_result.matched = [pw_result.matched;row,row;...
                            pw_result.annihilated(ind1,pw_result.annihilated(ind1,:)~=row) ...
                            pw_result.generated(ind2,pw_result.generated(ind2,:)~=row) ];
                        
                        pw_result.annihilated(ind1,:)=[];
                        pw_result.generated(ind2,:)=[];
                    end
                    alternative_matching(aii) = []; continue
                end
                
                % MG: Matched with a different pw and generated (happens
                %       because the location of the previous annihilation
                %       is stored, and new pws are matched to it even if a
                %       they still exist further down in the table)
                % [ 0 0  0 0 -P A ]       [ 0 0 0 0 0 A ]
                % [ P P -B 0 -A P ] -?->  [ P P 0 0 0 P ]
                % [-C B -P 0  0 0 ]       [ 0 0 0 0 0 0 ]
                % [-B C  C C  C C ]       [ 0 C C C C C ]
                
                if any(pw_result.matched(:,1)==row) && any(pw_result.generated(:)==row)
                    if correct_MG
                        ind = pw_result.matched(:,1)==row;
                        pw_result.new = [pw_result.new;pw_result.matched(ind,2)];
                        pw_result.matched(ind,:) = [];
                        
                        ind = any(pw_result.generated==row,2);
                        pw_result.new = [pw_result.new;pw_result.generated(ind,pw_result.generated(ind,:)~=row)];
                        pw_result.generated(ind,:)=[];
                        
                        pw_result.matched = [pw_result.matched; row,row];
                    end
                    alternative_matching(aii) = []; continue
                    
                elseif any(pw_result.matched(:,2)==row) && any(pw_result.annihilated(:)==row)
                    if correct_MG
                        ind = pw_result.matched(:,2)==row;
                        pw_result.lost = [pw_result.lost;pw_result.matched(ind,1)];
                        pw_result.matched(ind,:) = [];
                        
                        ind = any(pw_result.annihilated==row,2);
                        pw_result.lost = [pw_result.lost;pw_result.annihilated(ind,pw_result.annihilated(ind,:)~=row)];
                        pw_result.annihilated(ind,:)=[];
                        
                        pw_result.matched = [pw_result.matched; row,row];
                    end
                    alternative_matching(aii) = []; continue
                end
                
            end % end iteration over alternative matching
            
            % == Reconstruct the tracking_out table
            tracking_out = sortrows([...
                pw_result.matched;...
                pw_result.annihilated(:,1),-pw_result.annihilated(:,2);pw_result.annihilated(:,2),-pw_result.annihilated(:,1);...
                -pw_result.generated(:,1),pw_result.generated(:,2);-pw_result.generated(:,2),pw_result.generated(:,1);...
                pw_result.lost,-1i*pw_result.lost;...
                -1i*pw_result.new,pw_result.new...
                ]);
            
            %% PART 4: Check results and prepare output
            
            % Check that all pinwheels were added to the list using the
            % real IDs and not the row number
            if sum( histc( tracking(tracking_out(real(tracking_out(:,1))>0,1),1) , sort(tracking(real(tracking(:,1))>0,1)) )~=1 ) > 0
                error('The number of pinwheels changes when simplifying tracking array!')
            end
            if sum( histc( tracking(tracking_out(real(tracking_out(:,end))>0,end),end) , sort(tracking(real(tracking(:,end))>0,end)) )~=1 ) > 0
                error('The number of pinwheels changes when simplifying tracking array!')
            end
            
            % ===  Check errors (not in external, as no obj tracker is present)
            if obj.debug_mode > 0
                
                % Display remaining discrepancies in alternative matching
                for aii=1:length(alternative_matching)
                    [iy,ix]= find(tracking_out == alternative_matching(aii));
                    other_pw = abs(real(tracking_out(sub2ind(size(tracking_out),iy,mod(ix,2)+1))));
                    other_pw(other_pw==0)=[];
                    [iyy,~]= ind2sub(size(tracking_out),find(ismember(tracking_out(:),other_pw)));
                    original_matching = sortrows(real(tracking_out(unique([iy;iyy]),:)),[-1 -2]);
                    disp(['Alternative matching with row = ',num2str(alternative_matching(aii))])
                    disp(original_matching)
                end
                
                % If pinwheel x-y position available (i.e. internal tracking):
                if ~isempty(obj.debug_frames)
                    
                    % Distance of lost pinwheels to ROI boundary
                    idx= tracking(pw_result.new,end);
                    if ~isempty(idx)
                        dist = ceil(obj.debug_dist_ROI(sub2ind( size(obj.debug_dist_ROI), ceil(obj.debug_frames{end,1}.y(idx)) ,ceil(obj.debug_frames{end,1}.x(idx)) )))';
                        idx2=find(dist>obj.debug_limit);
                        if ~isempty(idx2)
                            disp(['Lost pinwheels ',mat2str(idx(idx2)),' are ',mat2str(dist(idx2)),' pixels away from ROI boundary.'])
                        end
                    end
                    
                    % Distance of new pinwheels to ROI boundary
                    idx= tracking(pw_result.lost,1);
                    if ~isempty(idx)
                        dist = ceil(obj.debug_dist_ROI(sub2ind( size(obj.debug_dist_ROI), ceil(obj.debug_frames{1,1}.y(idx)) ,ceil(obj.debug_frames{1,1}.x(idx)) )))';
                        idx2=find(dist>obj.debug_limit);
                        if ~isempty(idx2)
                            disp(['New pinwheels ',mat2str(idx(idx2)),' are ',mat2str(dist(idx2)),' pixels away from ROI boundary.'])
                        end
                    end
                end
            end
            
            %  ==  Convert row numbers to pinwheel ID
            if ~return_global_idx
                
                % --- loose LOST information (imaginary numbers)
                tracking_out = real(tracking_out);
                
                tmp = zeros(size(tracking_out));
                % matched pinwheels
                tmp(tracking_out(:,1)>0,1) = tmp(tracking_out(:,1)>0,1) + tracking(tracking_out(tracking_out(:,1)>0,1),1);
                tmp(tracking_out(:,1)<0,1) = tmp(tracking_out(:,1)<0,1) - tracking(abs(tracking_out(tracking_out(:,1)<0,1)),end);
                % annihilated pinwheels
                tmp(tracking_out(:,2)>0,2) = tmp(tracking_out(:,2)>0,2) + tracking(tracking_out(tracking_out(:,2)>0,2),end);
                tmp(tracking_out(:,2)<0,2) = tmp(tracking_out(:,2)<0,2) - tracking(abs(tracking_out(tracking_out(:,2)<0,2)),1);
                % pass
                tracking_out = tmp;
                
                % --- change labels in pw_result
                pw_result.matched = [real(tracking(pw_result.matched(:,1),1)) , real(tracking(pw_result.matched(:,2),end))];
                pw_result.lost = real(tracking(pw_result.lost,1));
                pw_result.new = real(tracking(pw_result.new,end));
                pw_result.annihilated = [real(tracking(pw_result.annihilated(:,1),1)) , real(tracking(pw_result.annihilated(:,2),1))];
                pw_result.generated = [real(tracking(pw_result.generated(:,1),end)) , real(tracking(pw_result.generated(:,2),end))];
                
            end
            
            % == reshape into output structure
            
            % First output: (for compatibility)
            % .ini all pinwheels that are matched from initial
            % .end the rest
            tmp = [];
            tmp.ini = tracking_out(tracking_out(:,1)>0,:);
            tracking_out(tracking_out(:,1)>0,:)=[];
            tmp.end = tracking_out(tracking_out(:,end)>0,:);
            tracking_out = tmp;
            
            if ~isempty(tracking_out.end)
                tracking_out.end = sortrows(tracking_out.end,2);
            end
            if ~isempty(tracking_out.ini)
                tracking_out.ini = sortrows(tracking_out.ini,1);
            end
            
        end
        
        %% ---------- PLOT MAPS FOR DEBUGGING
        function debug_gui(obj)
            % This GUI shows the maps in each interpolation step and the corresponding
            % pinwheels with its numbers. It helps to debug the code. Use the context
            % menu to change the mouse options:
            % Step: Scroll up/down to change interpolation step
            % Zoom: Scroll up/down to change zoom (don't open zoom menu with left click!)
            % Pan: Drag with left click to Pan the image
            % To finish the program close the figure.
            % chepe@nld.ds.mpg.de
            %%
            
            %% --------------------
            %      INITIATE GUI
            % ---------------------
            current_step = 1;
            text_is_row = false;
            
            % General square figure
            screenAspectRatio=get(0,'ScreenSize');
            screenAspectRatio=screenAspectRatio(4)/screenAspectRatio(3);
            
            % Table figure
            T.f = figure('Name','Tracking array','Units',...
                'Normalized','Position', [0.2 0.2 0.8*screenAspectRatio 0.8]);
            T.t = uitable('Units','Normalized','Parent', T.f, 'Position', [.05 .15 .9 .8]);
            set(T.t, 'Data', int16(obj.tracking_table));
            
            % Map figure
            H.f = figure('WindowScrollWheelFcn',@figScroll,...
                'Name','Pinwheel Tracking','Units',...
                'Normalized','Position', [0.1 0.1 0.9*screenAspectRatio 0.9],...
                'Color','w',...
                'menubar','none','numbertitle','off');
            
            % Map axes
            H.a = axes('Parent', H.f, 'Units','normalized','Position', [.05 .15 .9 .8],...
                'ButtonDownFcn',@buttonOption, ...
                'NextPlot', 'replacechildren',...
                'xlim',[1 size(obj.debug_frames{1,2},1)],'ylim',[0 size(obj.debug_frames{1,2},2)],...
                'xTick',[],'yTick',[],'Box','on','PlotBoxAspectRatio',[1 1 1]);
            
            % Interpolation slider
            H.sl = uicontrol('Style', 'slider','Units','normalized','Position',[0.25 0.09 0.5 0.05],...
                'Min',1,'Max',size(obj.debug_frames,1),'Value',1,'SliderStep',[1/size(obj.debug_frames,1) 10/size(obj.debug_frames,1)],'backgroundcolor',[1 1 1],...
                'Callback',{@sliderCallback,H});
            
            % Initial button
            H.b = uicontrol('Style', 'pushbutton', 'String', 'Initial',...
                'Units','Normalized','Position', [0.1 0.09 0.1 0.05],...
                'Callback',{@buttonIniCallback,H.sl});
            
            % End button
            H.b = uicontrol('Style', 'pushbutton', 'String', 'End',...
                'Units','Normalized','Position', [0.8 0.09 0.1 0.05],...
                'Callback', {@buttonEndCallback,H.sl});
            
            % show unique ID
            H.p = uicontrol('Style','checkbox', 'String', 'Show unique IDs',...
                'Units','Normalized','Position', [0.12 0.955 0.2 0.03],...
                'Callback',@checkID,'Backgroundcolor','w','Value',false);
            if isempty(obj.tracking_table)
                set(H.p,'enable','off')
            end
            
            % % context menu
            uimenu(H.f,'Label','Step','Callback', 'pan off;zoom off;');
            uimenu(H.f,'Label','Zoom','Callback', 'pan off;zoom on;');
            uimenu(H.f,'Label','Pan','Callback', 'zoom off;pan on;');
            
            % do initial plot
            doPlot
            
            % ---- Function that does the actual plotting
            function doPlot
                cla
                
                % colormap
                cm=[1                , 0.352941176470588, 0                 ;...
                    1                , 0.588235294117647, 0                 ;...
                    1                , 0.823529411764706, 0                 ;...
                    1                , 1                , 0                 ;...
                    0.823529411764706, 1                , 0                 ;...
                    0.549019607843137, 1                , 0                 ;...
                    0.156862745098039, 1                , 0                 ;...
                    0                , 1                , 0                 ;...
                    0                , 0.784313725490196, 0.196078431372549 ;...
                    0                , 0.549019607843137, 0.588235294117647 ;...
                    0                , 0.313725490196078, 1                 ;...
                    0                , 0                , 1                 ;...
                    0.392156862745098, 0                , 1                 ;...
                    0.705882352941177, 0                , 1                 ;...
                    1                , 0                , 0.478431372549020 ;...
                    1                , 0                , 0                 ];
                
                data = obj.debug_frames{current_step,2};
                ROI = ~isnan(data);
                
                ori=(angle(data))/(2*pi);ori(ori<0)=1+ori(ori<0);
                ori=ori(:);
                
                % find to which color interval it belongs
                intervals=linspace(0,1,2*size(cm,1)+1);
                [~,ori] = histc(ori,intervals(2:2:end));
                ori=ori+1;
                
                % match colors
                anglePlot = (cm(ori(:),:));
                anglePlot = reshape(anglePlot, [size(data,1),size(data,2),3]);
                
                fig2plot=repmat(double(ROI),[1 1 3]).*anglePlot;
                
                h=image(fig2plot);
                set(gca,'Ydir','reverse')
                hold on
                plot(obj.debug_frames{current_step,1}.x(obj.debug_frames{current_step,1}.sign==1),obj.debug_frames{current_step,1}.y(obj.debug_frames{current_step,1}.sign==1),'ko','MarkerFaceColor','w','MarkerSize',8)
                plot(obj.debug_frames{current_step,1}.x(obj.debug_frames{current_step,1}.sign==-1),obj.debug_frames{current_step,1}.y(obj.debug_frames{current_step,1}.sign==-1),'k^','MarkerFaceColor','w','MarkerSize',8)
                for ii=1:length(obj.debug_frames{current_step,1}.x)
                    x = obj.debug_frames{current_step,1}.x(ii);
                    y = obj.debug_frames{current_step,1}.y(ii);
                    if text_is_row
                        string = num2str(find(obj.tracking_table(:,current_step)==ii));
                    else
                        string = num2str(ii);
                    end
                    text(x+2,y+2,string,'Color','k','FontSize',16)
                end
                title(['Step ',num2str(current_step)],'FontSize',20)
                set(h,'HitTest','off') % mouse clicks go to axes category
                
            end
            
            % ---- Function that does the scrolling between steps
            function figScroll(src,evnt)
                
                if evnt.VerticalScrollCount > 0
                    %clc;disp('scroll down')
                    current_step=current_step+1;
                    if current_step>size(obj.debug_frames,1)
                        current_step=size(obj.debug_frames,1);
                    end
                    doPlot
                elseif evnt.VerticalScrollCount < 0
                    %clc;disp('scroll up')
                    current_step=current_step-1;
                    if current_step<1
                        current_step=1;
                    end
                    doPlot
                end
            end
            
            % ---- Functions for initial and final map buttons and interpolation slider
            function buttonIniCallback(h,~,H)
                current_step = 1;
                set(H,'value',current_step)
                doPlot
            end
            
            function buttonEndCallback(h,~,H)
                current_step = size(obj.debug_frames,1);
                set(H,'value',current_step)
                doPlot
            end
            
            function sliderCallback(varargin)
                [p,H] = varargin{[1,3]};  % calling handle and data structure.
                current_step=round(get(p,'value'));
                doPlot
            end
            
            % ------ Function to change which pinwheel ID is shown
            function checkID(h,~)
                if get(h,'Value')
                    text_is_row = true;
                else
                    text_is_row = false;
                end
                doPlot
            end
        end
        
    end
    
end

%% ==================  Class-Related Functions ==========================

%% ---------------------- STRUCTURE HELPER FUNCTIONS ---------------------
function s_new = structure_create(s,varargin)
% take given structure s as example and create new empty one with same fields
% Extra input parameters are new fields

% get field names
names = fieldnames(s);
% create new structure
s_new = struct;
% include fields
for ii=1:length(names)
    s_new.(names{ii})=[];
end
% add extra parameters
for ii=1:length(varargin)
    s_new.(varargin{ii})=[];
end

end

function s = structure_delete(s,indices)
% take structure s and delete the given indices in all fields

% get field names
names = fieldnames(s);
% delete indices in each field
for ii=1:length(names)
    s.(names{ii})(indices)=[];
end
end

function s = structure_extract(s,indices)
% take structure s and make a new structure with only the given indices in all fields

% get field names
names = fieldnames(s);
% extract indices in each field
for ii=1:length(names)
    s.(names{ii})=s.(names{ii})(indices);
end
end

function s = structure_join(s1,s2)
% Take two structures s1 and s2 and join in s combining the fields

% get field names
names = fieldnames(s1);
% combine in new structure
for ii=1:length(names)
    s.(names{ii})=[s1.(names{ii});s2.(names{ii})];
end
end

%% ----------- HELPER FUNCTION FOR PINWHEEL FINDER ------------------------
function [crossing,m,n] = line_line_intersection(P1x,P1y,P2x,P2y,Q1x,Q1y,Q2x,Q2y)
% function to calculate crossing between two lines (P and Q), defined each by
% two points (1 and 2) with its coordinates (x and y)
% http://en.wikipedia.org/wiki/Line-line_intersection

divisor=(P1x-P2x).*(Q1y-Q2y)-(P1y-P2y).*(Q1x-Q2x);
m=((P2y-P1y).*(Q1x-P1x)-(P2x-P1x).*(Q1y-P1y))./divisor;
n=((Q1x-Q2x).*(Q1y-P1y)-(Q1y-Q2y).*(Q1x-P1x))./divisor;

crossing=((0 <= m) & (m <= 1) & (0<=n) & (n<=1));

end

%% ------------ MATCH PINWHEELS HELPER FUNCTIONS --------------------------
function adjacency = test_adjacency(pw_distance,possible_adjacency,type,call_label)
% Function that takes an array of possible adjacency and generates an
% optimal adjacency model using the distance as cost. For this it separates
% the problem using 'separate_adjacency_array' into independent blocks and
% when the solution is not unique it uses 'construct_model' to find
% the combination that minimizes the total distance.
%
% INPUT:
% pw_distance = euclidian distance between pinwheels (cost function)
% possible_adjacency = logical array showing possible connections
% type = 'asymmentric' or 'symmetric' depending if the nodes in the columns
%       and rows are the same
% call_label = to display a warning when a sub_block is too large

% Start with empty adjacency matrix
adjacency = false(size(possible_adjacency));

% Separate problem into non-overlapping blocks
blocks = separate_adjacency_array(possible_adjacency,type);
for block_num=1:size(blocks,1)
    
    % If one of the sides of the block is empty, skip since there is no match
    if isempty(blocks{block_num,1}) || isempty(blocks{block_num,2})
        continue
    end
    
    % If only one model is possible in the block, add directly to adjacency
    if length(blocks{block_num,1})==1 && length(blocks{block_num,2})==1
        adjacency(blocks{block_num,1},blocks{block_num,2})=true;
        continue
    end
    
    % Create a model with maximal connectivity and minimal total distance
    pw_distance_block = pw_distance(blocks{block_num,1},blocks{block_num,2});
    possible_adjacency_block = possible_adjacency(blocks{block_num,1},blocks{block_num,2});
    
    if max(size(possible_adjacency_block))>15
        disp(['WARNING: large possible adjacency block N=',num2str(max(size(possible_adjacency_block))),' label=',call_label])
    end
    % If type is symmetric, use upper triangular to remove double labelling
    % while doing the model (distances are the same i->j and j->i )
    if strcmp(type,'symmetric')
        possible_adjacency_block = triu(possible_adjacency_block);
    end
    
    % find optimal adjacency model and include to adjacency array
    model = construct_model(possible_adjacency_block,pw_distance_block,type);
    
    % add model to adjacency matrix
    row_tmp = blocks{block_num,1}(model>0);
    col_tmp = blocks{block_num,2}(model(model>0));
    adjacency_idx=sub2ind(size(adjacency),row_tmp(:),col_tmp(:));
    adjacency(adjacency_idx)=true;
end

% Convert back to asymmetric (doubled matching label is used afterwards)
if strcmp(type,'symmetric')
    adjacency = adjacency+adjacency';
end
end

function blocks = separate_adjacency_array(possible_adjacency,type)
% Separate independent problems to check from the adjacency matrix. The
% code block-diagonalizes the array and returns it as independent blocks.
% This makes the process of trying the different optimal models much
% faster, as it only has to loop between the competing nodes
% INPUT:
% possible_adjacency= logical array of possible adjacent nodes
% type = 'asymmentric' or 'symmetric' depending if the nodes in the columns
%       and rows are the same
% OUTPUT:
% blocks = Cell giving the indices of the different separated blocks. The
%       number of rows corresponds to the number of blocks. In each row, the
%       first column has the row indices, the second column the column indices.

% Make a list of all columns that have to be sorted into blocks
cols_not_in_group = 1:size(possible_adjacency,2);

% Allocate initial block cell
blocks=cell(size(possible_adjacency,2),2);
block_ind=0;

while ~isempty(cols_not_in_group)
    block_ind=block_ind+1;
    
    
    % seed new block with first column that has not been sorted
    cols_to_check = cols_not_in_group(1);
    while ~isempty(cols_to_check)
        
        % find corresponding rows -> corresponding columns
        cols_to_check_next = [];
        % asymmetric case: rows and cols are not the same. Include
        % separately!
        if strcmp(type,'asymmetric')
            rows_checked = [];
            for col_ind=cols_to_check
                rows_to_check = find(possible_adjacency(:,col_ind))';
                for row_ind=rows_to_check
                    cols_to_check_next = unique([cols_to_check_next,find(possible_adjacency(row_ind,:))]);
                end
                rows_checked = [rows_checked,rows_to_check];
            end
            
            % symmetric case: rows and cols are the same. Every included row is
            % added as a column and vice versa
        elseif strcmp(type,'symmetric')
            for col_ind=cols_to_check
                rows_to_check = [cols_to_check,find(possible_adjacency(:,col_ind))'];
                for row_ind=rows_to_check
                    cols_to_check_next = unique([cols_to_check_next,find(possible_adjacency(row_ind,:))]);
                end
            end
            rows_checked = cols_to_check_next;
        end
        
        % add checkd indices in current block
        blocks{block_ind,1} = union(blocks{block_ind,1},rows_checked);
        blocks{block_ind,2} = union(blocks{block_ind,2},cols_to_check);
        
        % remove indices that already make part of the block
        cols_to_check = setdiff(cols_to_check_next,blocks{block_ind,2});
    end
    
    % remove added columns from the unsorted list
    cols_not_in_group = setdiff(cols_not_in_group,blocks{block_ind,2});
end

% clean by deleting the rest of group list
blocks(block_ind+1:end,:)=[];

end

function model = construct_model(possible_adjacency,costs_array,type)
% Function that finds the optimal adjacency model for the given possible
% connections using a given cost_array. The function starts a model by scanning
% through the rows of the possible_adjacency array, picks a column, 'zeros'
% this column/row to make a new array and sends it RECURSIVELY to a copy of
% the same function to complete the model. Once the new model is filled the NESTED function
% 'compare_model_cost' calculates the cost per connection and if it is lower than the cost
% of the saved model, exchanges it by the new model.
%
% INPUT:
% possible_adjacency = logical array of possible adjacent nodes
% costs_array = cost of the possible connections (e.g. euclidean distance)
% type = if the nodes are 'symmetric' (the same nodes are in the rows
%           and the columns) or 'asymmetric'
%
% OUTPUT:
% model = vector representing the adjacency array. The indices correspond to
%       the rows of the array and the value inside to the column with 1=adjacent.
%       If the index is filled with 0, then there is no adjacent node.

% add one column to 'costs_array' to calculate cost since the model can have zeros, i.e no connection
costs_expanded=[zeros(size(costs_array,1),1),costs_array];

% initiate adjacency model and cost
model = zeros(size(possible_adjacency,1),1);
cost = Inf;

% initiate new model
add_row_to_model(model,possible_adjacency)

% -- NESTED FUNCTION to add one line to the model. If empty, calculate
%       final cost of the model
    function add_row_to_model(model_current,adjacency_current)
        % if only zeros in adjacency matrix, calculate cost of model and return
        if sum(adjacency_current(:)==1)==0
            compare_model_cost(model_current)
            return
        else
            % check possible models in a recursive loop (go only through possibly adjacent)
            ind_row = find(sum(adjacency_current,2)>0,1,'first');
            % find filled columns in the identified row
            for ind_col = [find(adjacency_current(ind_row,:)>0) 0] % add a zero to have no match in row
                % add row to model
                model_next = model_current;
                model_next(ind_row)=ind_col;
                % reduce adjacency
                adjacency_next = adjacency_current;
                adjacency_next(ind_row,:)=0;
                if ind_col>0
                    adjacency_next(:,ind_col)=0;
                end
                if strcmp(type,'symmetric')
                    % remove also corresponding row and column
                    if ind_col>0
                        adjacency_next(ind_col,:)=0;
                    end
                    adjacency_next(:,ind_row)=0;
                end
                % calculate next row in the model
                add_row_to_model(model_next,adjacency_next);
                
                % if symmetric, do also col=row
                if strcmp(type,'symmetric') && ind_col>0
                    ind_row_s = ind_col;
                    for ind_col_s = find(adjacency_current(ind_row_s,:)>0)
                        % add row to model
                        model_next = model_current;
                        model_next(ind_row_s)=ind_col_s;
                        % reduce adjacency
                        adjacency_next = adjacency_current;
                        adjacency_next(ind_row_s,:)=0; adjacency_next(:,ind_row_s)=0;
                        adjacency_next(:,ind_col_s)=0; adjacency_next(ind_col_s,:)=0;
                        % calculate next row in the model
                        add_row_to_model(model_next,adjacency_next);
                    end
                end
            end
        end
    end

% -- NESTED FUNCTION to calculate cost per connection of input model and
%       if smaller than current cost, use it as base model
    function compare_model_cost(model_to_compare)
        % only if model is not empty
        if sum(model_to_compare>0)>=sum(model>0) % only use models that match more pinwheels
            cost_to_compare = ...
                sum(costs_expanded(sub2ind(size(costs_expanded), [1:length(model)]', model_to_compare+1)))/sum(model_to_compare>0)^2;
            if cost_to_compare<cost
                model = model_to_compare;
                cost = cost_to_compare;
            end
        end
    end
end

%% ------------- TRASH FUNCTIONS -----------------------------------------
% function indices = structure_find(s,varargin)
% % Find indices in structure s where the pairs ['fieldname',value] in varargin are found
%
% % assume all indices satisfy the test (take the first test argument to find size)
% indices = true(size(s.(varargin{1})));
% % remove indices applying the tests given in varargin
% for test = 1:length(varargin)/2
%     indices(~ismember(s.(varargin{2*test-1}),varargin{2*test}))=false;
% end
% end
% function [possible_adjacency,adjacency]=reduce_adjacency_array_2(possible_adjacency)
% % Check pairs that are unambiguous (only one 1 in row x column pair), remove
% % from possible_adjacency and add to adjacency array.
% % INPUT:
% % possible_adjacency: logical array of possible adjacent nodes
% % type = 'asymmentric' or 'symmetric' depending if the nodes in the columns
% %       and rows are the same
% % OUTPUT:
% % possible_adjacency: reduced array where unambiguous pairs have been removed
% % adjacency: logical adjacency array
%
% % initiate adjacency array
% adjacency = false(size(possible_adjacency));
%
% % find rows with only one 1
% for row_id=find(sum(possible_adjacency,2)==1)'
%     % find the column of the 1
%     col_id = find(possible_adjacency(row_id,:));
%     % if that column also only has one 1, then the pair is unique
%     if sum(possible_adjacency(:,col_id))==1
%         adjacency(row_id,col_id)=true;
%         possible_adjacency(row_id,:)=false;
%         possible_adjacency(:,col_id)=false;
%     end
% end
%
% end
% function model = construct_optimal_model(possible_adjacency,costs_array,type)
% % Function that finds the optimal adjacency model for the given possible
% % connections using a given cost_array. The function starts a model by scanning
% % through the rows and columns of the possible_adjacency array, removes
% % this column/row to make a new array and sends it RECURSIVELY to a copy of
% % the same function to complete the model. Once the new model is filled (the size
% % of the model depends on the level of recursion) the NESTED function
% % 'compare_model_cost' calculates the cost per connection and if it is lower than the cost
% % of the saved model, exchanges it by the new model.
% %
% % INPUT:
% % possible_adjacency = logical array of possible adjacent nodes
% % costs_array = cost of the possible connections (e.g. euclidean distance)
% % type = if the nodes are 'symmetric' (the same nodes are in the rows
% %           and the columns) or 'asymmetric'
% %
% % OUTPUT:
% % model = vector representing the adjacency array. The indices correspond to
% %       the rows of the array and the value inside to the column with 1=adjacent.
% %       If the index is filled with 0, then there is no adjacent node.
%
% % add one column to 'costs_array' to calculate cost since the model can have zeros, i.e no connection
% costs_expanded=[zeros(size(costs_array,1),1),costs_array];
%
% % get size of input possible_adjacency
% [rows,cols] = size(possible_adjacency);
%
% % initiate adjacency model and cost
% model = zeros(rows,1);
% cost = Inf;
%
% % check possible models in a recursive loop (go only through possibly adjacent)
% for ind_row = find(sum(possible_adjacency,2)>0)'
%
%     % find filled columns in the identified row
%     for ind_col = find(possible_adjacency(ind_row,:)>0)
%
%         % START NEW MODEL and add the starting index
%         model_current = zeros(rows,1);
%         model_current(ind_row)=ind_col;
%
%         % fix the values in ind_row and ind_col and reduce the possible_adjacency array
%         reduced_rows = 1:rows;
%         reduced_cols = 1:cols;
%         switch type
%             case 'asymmetric'
%                 reduced_rows(ind_row)=[];
%                 reduced_cols(ind_col)=[];
%             case 'symmetric' % remove also corresponding row
%                 reduced_rows([ind_row ind_col])=[];
%                 reduced_cols([ind_row ind_col])=[];
%         end
%         reduced_possible_adjacency = possible_adjacency(reduced_rows,reduced_cols);
%         reduced_cost_array = costs_array(reduced_rows,reduced_cols);
%
%
%         % check if the reduced model needs to be recomposed recursively to
%         % find the one with lowest cost or if it only defines an unique
%         % model
%         if ~isempty(reduced_possible_adjacency)
%
%             % make from symmetric -> asymetric to only count matches in rows
%             switch type
%                 case 'asymmetric'
%                     test = reduced_possible_adjacency;
%                 case 'symmetric'
%                     test= reduced_possible_adjacency+reduced_possible_adjacency';
%             end
%             test = sum(test,2);
%
%             if sum(test>1)==0  % only 1 to 1 mapping
%                 % --- 1) add matching pairs to model
%                 [row_match,col_match]=find(reduced_possible_adjacency==1);
%                 model_current(reduced_rows(row_match))=col_match;
%
%             else % double matching -> reduce model recursively
%
%                 %-- 1) Call the function with the reduced array !RECURSIVE CALL OF FUNCTION!
%                 model_recieved = construct_optimal_model(reduced_possible_adjacency,reduced_cost_array,type);
%
%                 %-- 2) Complement the current model with the recieved model
%                 % shift the values in the recieved model to accomodate deleted column(s) when reducing
%                 switch type
%                     case 'asymmetric'
%                         model_recieved(model_recieved>=ind_col)=model_recieved(model_recieved>=ind_col)+1;
%                     case 'symmetric'
%                         cols2add = sort([ind_col,ind_row],'ascend');
%                         model_recieved(model_recieved>=cols2add(1))=model_recieved(model_recieved>=cols2add(1))+1;
%                         model_recieved(model_recieved>=cols2add(2))=model_recieved(model_recieved>=cols2add(2))+1;
%                 end
%                 % add the recieved model to the current model
%                 model_current(reduced_rows) = model_recieved;
%             end
%         end
%
%         % Check the model now that it is complete !LOWEST LEVEL OF RECURSION!
%         compare_model_cost(model_current)
%     end
% end
%
% % -- NESTED FUNCTION to calculate cost per connection of input model and if smaller than current cost, use it as base model
%     function compare_model_cost(model_to_compare)
%         % only if model is not empty
%         if sum(model_to_compare)>0
%             cost_to_compare = ...
%                 sum(costs_expanded(sub2ind(size(costs_expanded), [1:length(model)]', model_to_compare+1)))/sum(model_to_compare>0)^2;
%             if cost_to_compare<cost
%                 model = model_to_compare;
%                 cost = cost_to_compare;
%             end
%         end
%     end
% end

%% OLD REDUCE TRACKING

% function tracking =reduce_tracking(obj,tracking)
% % This code reads a tracking array and makes a reduced version of it showing
% % which pinwheels correspond to which (or gets annihilated/generated lost/new)
% % There are 5 possible cases. The pinwheel can:
% % 1) be matched from beginning to end
% % 2) get lost during the interpolations
% % 3) get annihilated
% % 4) appear during the interpolations
% % 5) get generated
% % For the cases 3) and 5), sometimes a pinwheel gets annihilated/generated
% % by a member of a pair that survives until the end/beginning of the
% % interpolation. In this case that pinwheel is tracked. The
% % generation/annihilation of pairs can also happen multiple times during
% % the interpolations, such that this check has to be done recursively
% 
% if nargin==1
%     tracking = double(obj.tracking_table);
% end
% %% % Clean tracking by removing spurious new->lost pinwheels (i.e. pinwheels
% % that appear and disappear without a partner)
% 
% tracking(:,size(tracking,2)+1) = NaN;
% tracking = tracking.';
% for spurious_length = 1: (size(tracking,2)-2)
%     idx =   strfind(sign(real(tracking(:))+imag(tracking(:))).',[0 ones(1,spurious_length) 0]);
%     for toDelete=1:length(idx)
%         tracking(idx(toDelete)+[0:spurious_length])=0;
%     end
% end
% tracking = tracking.';
% tracking(:,end) = [];
% % tracking = int16(tracking);
% 
% %% Start empty array and go row by row
% tracking_out = [];
% try
%     for row=1:size(tracking,1)
%         
%         % CASE 0: Phatonm pinwheels appearing only during interpolation steps
%         if real(tracking(row,1))<1 && real(tracking(row,end))<1
%             continue
%         end
%         
%         % CASE 1: Pinwheel is matched from beginning to end
%         if real(tracking(row,1))>0 && real(tracking(row,end))>0
%             tracking_out = [tracking_out;real(tracking(row,1)) real(tracking(row,end))];
%             continue
%         end
%         
%         % CASE 2: Pinwheel gets lost at some point in the interpolation
%         % [control for annihilated->re-generated pinwheels]
%         if real(tracking(row,1))>0 && real(tracking(row,end))==0 && ...
%                 sum(tracking(row,find(real(tracking(row,:))>0,1,'last')+1:end))==0
%             tracking_out = [tracking_out;real(tracking(row,1)) 0];
%             continue
%         end
%         
%         % CASE 3: Pinwheel gets annihilated
%         if real(tracking(row,1))>0 && real(tracking(row,end))<1 && ...
%                 find(imag(tracking(row,:))<0,1,'last') > find(real(tracking(row,:))>0,1,'last')
%             
%             % initiate annihilation - generation loop
%             row_loop = row;
%             continue_loop = true;
%             while continue_loop
%                 
%                 % find co-annihilated
%                 test_when = find(imag(tracking(row_loop,:))<0,1,'last');
%                 test_with = find(imag(tracking(:,test_when))==-real(tracking(row_loop,test_when-1)) & real(tracking(:,test_when-1))>0);
%                 
%                 % check if the co-annihilated partner exists from the beginning
%                 if real(tracking(test_with,1))>0
%                     tracking_out = [tracking_out;real(tracking(row,1)) -real(tracking(test_with,1))];
%                     continue_loop = false;
%                     
%                     % check if the co-annihilated partner is generated
%                 elseif find(real(tracking(test_with,:))<0,1,'first') < find(real(tracking(test_with,:))>0,1,'first')
%                     
%                     % find co-generated
%                     test_when = find(real(tracking(test_with,:))<0,1,'first');
%                     test_with = find(real(tracking(:,test_when))==-real(tracking(test_with,test_when+1)) & real(tracking(:,test_when+1))>0);
%                     
%                     % check if the co-generated partner exists until the end
%                     if real(tracking(test_with,end))>0
%                         tracking_out = [tracking_out;real(tracking(row,1)) real(tracking(test_with,end))];
%                         continue_loop = false;
%                         
%                         % check if co-generated partner is annihilated
%                     elseif find(imag(tracking(test_with,:))<0,1,'last') > find(real(tracking(test_with,:))>0,1,'last')
%                         row_loop = test_with;
%                         continue_loop = true;
%                         
%                         % otherwise matching pinwheel is lost
%                     else
%                         tracking_out = [tracking_out;real(tracking(row,1)) -real(tracking(row,1))];
%                         continue_loop = false;
%                     end
%                     
%                     % the matching pinwheel appears at some point
%                 else
%                     tracking_out = [tracking_out;real(tracking(row,1)) -real(tracking(row,1))];
%                     continue_loop = false;
%                 end
%             end
%             continue
%         end
%         
%         % CASE 4: Pinwheel appears at some point in the tracking
%         % [control for annihilated->re-generated pinwheels]
%         if real(tracking(row,1))==0 && real(tracking(row,end))>0 && ...
%                 sum(tracking(row,1:find(real(tracking(row,:))>0,1,'first')-1))==0
%             tracking_out = [tracking_out;0 real(tracking(row,end))];
%             continue
%         end
%         
%         % CASE 5: Pinwheel gets generated
%         if real(tracking(row,1))<1 && real(tracking(row,end))>0 && ...
%                 find(real(tracking(row,:))<0,1,'first') < find(real(tracking(row,:))>0,1,'first')
%             
%             % initiate generation - annihilation loop
%             row_loop = row;
%             continue_loop = true;
%             while continue_loop
%                 
%                 % find co-generated
%                 test_when = find(real(tracking(row_loop,:))<0,1,'first');
%                 test_with = find(real(tracking(:,test_when))==-real(tracking(row_loop,test_when+1)) & real(tracking(:,test_when+1))>0);
%                 
%                 % check if the co-generated partner exists from until the end
%                 if real(tracking(test_with,end))>0
%                     tracking_out = [tracking_out;-real(tracking(test_with,end)) real(tracking(row,end))];
%                     continue_loop = false;
%                     
%                     % check if the co-generated partner is annihilated
%                 elseif find(imag(tracking(test_with,:))<0,1,'last') > find(real(tracking(test_with,:))>0,1,'last')
%                     
%                     % find co-annihilated
%                     test_when = find(imag(tracking(test_with,:))<0,1,'last');
%                     test_with = find(imag(tracking(:,test_when))==-real(tracking(test_with,test_when-1)) & real(tracking(:,test_when-1))>0);
%                     
%                     % check if the co-annihilated partner exists from the
%                     % beginning
%                     if real(tracking(test_with,1))>0
%                         tracking_out = [tracking_out;real(tracking(test_with,1)) real(tracking(row,end))];
%                         continue_loop = false;
%                         
%                         % check if co-annihilated partner is generated
%                     elseif find(real(tracking(test_with,:))<0,1,'first') < find(real(tracking(test_with,:))>0,1,'first')
%                         row_loop = test_with;
%                         continue_loop = true;
%                         
%                         % otherwise matching pinwheel is lost
%                     else
%                         tracking_out = [tracking_out;-real(tracking(row,end)) real(tracking(row,end))];
%                         continue_loop = false;
%                     end
%                     
%                     % the partner appears at some point
%                 else
%                     tracking_out = [tracking_out;-real(tracking(row,end)) real(tracking(row,end))];
%                     continue_loop = false;
%                 end
%             end
%             continue
%         end
%     end
% catch err
%     % problem when reducing the tracking array
%     disp(['Unknown error when simplifying tracking array!'])
%     rethrow(err)
% end
% 
% 
% %% Double tracking is possible in generation/annihilation loops. Find unique
% tracking_out = unique(tracking_out,'rows');
% 
% % Check that all pinwheels were added to the list
% if (length(tracking_out(tracking_out(:,1)>0,1))~=length(tracking(real(tracking(:,1))>0,1))) ...
%         || (length(unique(tracking_out(tracking_out(:,1)>0,1)))~=length(tracking(real(tracking(:,1))>0,1)))
%     error(['The number of pinwheels changes when simplifying tracking array!'])
% end
% if (length(tracking_out(tracking_out(:,end)>0,end))~=length(tracking(real(tracking(:,end))>0,end))) ...
%         || (length(unique(tracking_out(tracking_out(:,end)>0,end)))~=length(tracking(real(tracking(:,end))>0,end)))
%     error(['The number of pinwheels changes when simplifying tracking array!'])
% end
% 
% %% Check errors (not in external, as no obj tracker is present)
% if obj.debug_mode > 0
%     
%     % Lost pinwheels
%     idx= tracking_out(tracking_out(:,2)==0 | (tracking_out(:,2)<0 & tracking_out(:,1)==-tracking_out(:,2)),1);
%     if ~isempty(idx)
%         dist = ceil(obj.debug_dist_ROI(sub2ind( size(obj.debug_dist_ROI), ceil(obj.debug_frames{1,1}.y(idx)) ,ceil(obj.debug_frames{1,1}.x(idx)) )))';
%         idx2=find(dist>obj.debug_limit);
%         if ~isempty(idx2)
%             disp(['Interpolating from first to last step the lost pinwheels ',mat2str(idx(idx2)),' are ',mat2str(dist(idx2)),' pixels away from ROI boundary.'])
%         end
%     end
%     
%     % New pinwheels
%     idx= tracking_out(tracking_out(:,1)==0 | (tracking_out(:,1)<0 & tracking_out(:,2)==-tracking_out(:,1)),2);
%     if ~isempty(idx)
%         dist = ceil(obj.debug_dist_ROI(sub2ind( size(obj.debug_dist_ROI), ceil(obj.debug_frames{end,1}.y(idx)) ,ceil(obj.debug_frames{end,1}.x(idx)) )))';
%         idx2=find(dist>obj.debug_limit);
%         if ~isempty(idx2)
%             disp(['Interpolating from first to last step the new pinwheels ',mat2str(idx(idx2)),' are ',mat2str(dist(idx2)),' pixels away from ROI boundary.'])
%         end
%     end
% end
% 
% %% reshape into output structure
% clear tracking
% tracking.ini = tracking_out(tracking_out(:,1)>0,:);
% tracking_out(tracking_out(:,1)>0,:)=[];
% tracking.end = tracking_out(tracking_out(:,end)>0,:);
% 
% if ~isempty(tracking.end)
%     tracking.end = sortrows(tracking.end,2);
% end
% if ~isempty(tracking.ini)
%     tracking.ini = sortrows(tracking.ini,1);
% end
% end
