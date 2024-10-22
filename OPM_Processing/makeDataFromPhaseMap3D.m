function data = makeDataFromPhaseMap3D(PhaseMap,stim_order)


    %% make complex
    if isreal(PhaseMap)
        z = exp(1i*PhaseMap);
    else
        z = PhaseMap;
    end

    trials = size(PhaseMap,1);
    MapDim = size(PhaseMap,[2 3]);
    data = zeros([ MapDim length(stim_order) trials]);

    for ii = 1:length(stim_order)
        stim = stim_order(ii);
        if ~isnan(stim)            
            for jj = 1: trials
                A =abs(z(jj,:,:)).*cos(angle(z(jj,:,:))-(2*pi*stim/180));
                data(:,:,ii,jj) = reshape(A,MapDim);
            end
        end
    end
end



% 
%     %% fix dimensions
%     if size(PhaseMap,1)<size(PhaseMap,3)
%         PhaseMap = reshape(PhaseMap,[size(PhaseMap,2) size(PhaseMap,3) size(PhaseMap,1)]);
%     end
% 
%     %% make complex
%     if isreal(PhaseMap)
%         z = exp(1i*PhaseMap);
%     else
%         z = PhaseMap;
%     end
% 
%     trials = size(PhaseMap,3);
%     data = zeros([size(PhaseMap,1) size(PhaseMap,2) length(stim_order) trials]);
% 
%     for ii = 1:length(stim_order)
%         stim = stim_order(ii);
%         if ~isnan(stim)            
%             for jj = 1: trials
%                 A =abs(z(:,:,jj)).*cos(angle(z(:,:,jj))-(2*pi*stim/180));
%                 data(:,:,ii,jj) = A;
%             end
%         end
%     end