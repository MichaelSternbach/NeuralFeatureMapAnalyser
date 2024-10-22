function Samples= getSamples(PwInfos,PwDensityType,data_obj,local)
    if local
        Samples = zeros(sum(data_obj.ROI(:)),length(PwInfos));
    else
        Samples = zeros(1,length(PwInfos));
    end
    
    for ii=1:length(PwInfos)
        
        if local
            Samples(:,ii) = data_obj.array2vector(PwInfos{ii}.(PwDensityType));
        else
            Samples(ii) = PwInfos{ii}.(PwDensityType);
        end
    end  
end