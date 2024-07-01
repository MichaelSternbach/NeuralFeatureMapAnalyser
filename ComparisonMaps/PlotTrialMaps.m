trials = 10:15;

z = make_map(mean(data_obj.data(:,:,:,trials),4),data_obj.info.stim_order,data_obj.ROI,false);
z_filtered =data_obj.filter_map(z,true);

figure
plot_map(z_filtered,data_obj.ROI,0,1)
title('Fermi filtered map ')

figure
selMap = abs(z_filtered.*data_obj.ROI);

colormap gray
imagesc(selMap/max(selMap,[],'all'))
title('selectivity map ')