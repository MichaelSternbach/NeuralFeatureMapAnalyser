function data_obj_rand = getRandomisedDataHandle(data_obj)
    data = randomizeData(data_obj.data);
    data_obj_rand =  data_handle_corrected(data_obj.info,data,data_obj.ROI);
end

