function data = generateModelData(data_info,trials,column_spacing_mm)
    map_dim = max(data_info.field_size_pix);
    [X,Y] = meshgrid(1:map_dim,1:map_dim);
    column_spacing_pix = column_spacing_mm*data_info.pix_per_mm;
    num_hc = map_dim/column_spacing_pix;
    [z,~] = makeOPM(data_info.animal,X,Y,num_hc);
    
    z = z(1:data_info.field_size_pix(1),1:data_info.field_size_pix(2));

    data = makeData(z,data_info.stim_order,trials);

end

function data = makeData(z,stim_order,trials)
    data = zeros([size(z) length(stim_order) trials]);

    for ii = 1:length(stim_order)
        stim = stim_order(ii);
        if ~isnan(stim)
            % A =abs(z.*conj(ones(size(z)).*exp(2j*pi*stim/180)));
            A =abs(z).*cos(angle(z)-(2*pi*stim/180));
            for jj = 1: trials
                data(:,:,ii,jj) = A;
            end
        end
    end
end