function PwDensityCI=getPwDensityCI(PwInfosBS,PwInfosJS,data_obj,PwDensityType,local,alpha)
    
    %% get PwDensits of jackknife samples
    SamplesBS= getSamples(PwInfosBS,PwDensityType,data_obj,local);

    %% get PwDensits of jackknife samples
    SamplesJS = getSamples(PwInfosJS,PwDensityType,data_obj,local);

    %% Calc CIs
    if local
        PwDensityCI_Vector = bootstrap_ci(SamplesBS,data_obj.array2vector(PwInfosBS{1}.(PwDensityType)),SamplesJS,alpha);
        PwDensityCI = zeros([size(data_obj.ROI) 3]);
        PwDensityCI(:,:,1) = data_obj.vector2array(PwDensityCI_Vector(:,1));
        PwDensityCI(:,:,2)=   PwInfosBS{1}.(PwDensityType);
        PwDensityCI(:,:,3) = data_obj.vector2array(PwDensityCI_Vector(:,2));
    else
        PwDensityCI = bootstrap_ci(SamplesBS,PwInfosBS{1}.(PwDensityType),SamplesJS,alpha);
    end

end

% function Samples= getSamples(PwInfos,PwDensityType,data_obj,local)
%     if local
%         Samples = zeros(sum(data_obj.ROI(:)),length(PwInfos));
%     else
%         Samples = zeros(1,length(PwInfos));
%     end
%     
%     for ii=1:length(PwInfos)
%         
%         if local
%             Samples(:,ii) = data_obj.array2vector(PwInfos{ii}.(PwDensityType));
%         else
%             Samples(ii) = PwInfos{ii}.(PwDensityType);
%         end
%     end  
% end