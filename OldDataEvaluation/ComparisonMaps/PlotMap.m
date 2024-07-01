z_filtered =data_obj.filter_map(data_obj.read_map(80),true);

figure
plot_map(z_filtered)
title('Fermi filtered map')