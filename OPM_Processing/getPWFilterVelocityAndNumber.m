function [PwFilterVelocity,PwFilterVelocityMean,PWNumber,lowpass_cutoffs] = getPWFilterVelocityAndNumber(data_obj,tracker_obj,min_lowpass_mm,max_lowpass_mm,n_steps_lowpass)
    orig_lowpass = data_obj.filter_parameters.lowpass;
    
    %% define lowpass_cutoffs
    
    lowpass_cutoffs = linspace(min_lowpass_mm,max_lowpass_mm,n_steps_lowpass);
    
    
    %% prep variables
    
    PwFilterVelocityMean = zeros([1,n_steps_lowpass-1]);
    PWNumber = zeros([1,n_steps_lowpass]);
    
    %% make 1 maps
    data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(1))
    z_1 = data_obj.filter_map(data_obj.read_map());
    pinwheels_1 =tracker_obj.find_pinwheels(z_1,data_obj.ROI);
    PWNumber(1) = getPWNumer(pinwheels_1);
    
    %% iterate through cut offs
    
    for ii = 2:n_steps_lowpass
        
        %% make maps
        data_obj.set_filter_parameters('lowpass',lowpass_cutoffs(ii))
        z_2 = data_obj.filter_map(data_obj.read_map());
        pinwheels_2 =tracker_obj.find_pinwheels(z_2,data_obj.ROI);
        PWNumber(ii) = getPWNumer(pinwheels_2);
        
        %% track pinwheels
        tracking = tracker_obj.interpolate(z_1,z_2,data_obj.ROI);
        
        %% calclate PW velocity
        PwFilterVelocity{ii-1} = getPwPositionDifference(tracking,pinwheels_1,pinwheels_2);%/abs(lowpass_cutoffs(ii)-lowpass_cutoffs(ii-1));
        PwFilterVelocityMean(ii-1) = mean(PwFilterVelocity{ii-1},'all');
        
        z_1 = z_2;
        pinwheels_1 =pinwheels_2;
        

    end
    
    %% revert to original lowpass
    data_obj.set_filter_parameters('lowpass',orig_lowpass)
end


function PWDistances = getPwPositionDifference(tracking,pinwheels_1,pinwheels_2)
    PWNumber = getPWNumer(pinwheels_1);
    PWDistances = [];
    for ii = 1:PWNumber
       pw1 = tracking.ini(ii,1);
       x1 =  pinwheels_1.x(pw1);
       y1 =  pinwheels_1.y(pw1);
       
       pw2 = tracking.ini(pw1,2);
       if pw2 ~= 0
           x2 =  pinwheels_2.x(pw2);
           y2 =  pinwheels_2.y(pw2);

           PWDistances = [PWDistances,sqrt((x1-x2)^2+(y1-y2)^2)];
       end
    end
    %MeanPwPositionDifference = mean(PWDistances,'all');
end

function PWNumber = getPWNumer(pinwheels)
    PWNumber = size(pinwheels.label,1);
end