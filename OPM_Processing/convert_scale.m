function scale = convert_scale(scale,data_obj)
    if scale > 1
        size_map = max(data_obj.info.field_size_pix);
        if scale > size_map
            scale = 1;
        else
            scale = scale/size_map;
        end
    end
end