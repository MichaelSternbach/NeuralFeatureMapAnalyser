function z = makeOPMFromPhaseMap3D(PhaseMap)


    %% make complex
    if isreal(PhaseMap)
        z = exp(1i*PhaseMap);
    else
        z = PhaseMap;
    end

    z = transform_array_format(z);

    

%     trials = size(PhaseMap,2);
%     MapDim = size(PhaseMap,[1 3 4]);
%     data = zeros([ MapDim length(stim_order) trials]);
% 
%     for ii = 1:length(stim_order)
%         stim = stim_order(ii);
%         if ~isnan(stim)            
%             for jj = 1: trials
%                 A =abs(z(jj,:,:)).*cos(angle(z(jj,:,:))-(2*pi*stim/180));
%                 data(:,:,ii,jj) = reshape(A,MapDim);
%             end
%         end
%     end
end

function transformed_array = transform_array_format(input_array)
    % Function to transform a 4D array from format [Z, T, X, Y] to [Z, X, Y, T]
    %
    % input_array: 4D array with dimensions [Z, T, X, Y]
    % Output:
    % transformed_array: 4D array with dimensions [Z, X, Y, T]
    
    % Rearrange the dimensions using permute
    transformed_array = permute(input_array, [1, 3, 4, 2]);
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