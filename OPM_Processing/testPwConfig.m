function [selectivities_pw_all,selectivities_pw_rand_all] = testPwConfig(data_obj,bootstrapsamples,ResultDataFolder)
    
    %% get pinwheel stats
    
    DataFile = [ResultDataFolder 'PwConfig.mat'];
    if ~isfile(DataFile)

        
        %% prepare randomized data
        data_rand = randomizeData(data_obj.data);
        data_obj_rand =  data_handle_corrected(data_obj.info,data_rand,data_obj.ROI);
        if data_obj.GIF_apply
            data_obj_rand.activateGIF(true,data_obj.SN_TH)
        end
        if data_obj.lsm_applied
            data_obj_rand.apply_LSM(true)
        end
        data_obj_rand.prepare_samples_array(bootstrapsamples)

        %% find pinwheels mean map
        z = data_obj.filter_map(data_obj.read_map);
        ROI = data_obj.ROI;
        [~,~,~,pw_pos.x,pw_pos.y,~,~] = find_pinwheels(z,0,ROI);
    

        %% find sum selectivities pinwheels
        selectivities_pw_all = SumSelectivitiesPw(pw_pos,data_obj);
        selectivities_pw_rand_all = SumSelectivitiesPw(pw_pos,data_obj_rand);
    
        save(DataFile,'data_obj_rand', ...
            "selectivities_pw_all",'selectivities_pw_rand_all')
    else
        load(DataFile,'data_obj_rand',...
            "selectivities_pw_all",'selectivities_pw_rand_all')
    end


    %% normalize selectivities by selectivity of mean map
    z = data_obj.filter_map(data_obj.read_map);
    z_rand = data_obj_rand.filter_map(data_obj_rand.read_map);
    ROI = data_obj.ROI;

    selectivities_pw_all = 1-selectivities_pw_all./nanmean(abs(z(ROI)));
    selectivities_pw_rand_all = 1-selectivities_pw_rand_all./nanmean(abs(z_rand(ROI)));

end

function selectivities_pw = SumSelectivitiesPw(pw_pos,data_obj)
    
    selectivities_pw = zeros([1 size(data_obj.samples_array,3)]);
    for ii = 2:size(data_obj.samples_array,3)
        S = 0;
        z = data_obj.filter_map(data_obj.read_map(ii));
        %z = (z-mean(z,'all'))./std(z);
        for jj = 1:size(pw_pos.x,1)
            x = round(pw_pos.x(jj));
            y = round(pw_pos.y(jj));
            S = S + abs(z(y,x));
        end
        selectivities_pw(ii) = S/size(pw_pos.x,1);
    end
    %selectivities_pw = selectivities_pw(2:end);
end