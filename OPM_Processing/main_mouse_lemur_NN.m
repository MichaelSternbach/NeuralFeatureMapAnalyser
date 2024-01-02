load microcebus_Huber.mat

plot_nn_dist('d', d, 'microcebus_nn_dist_d')
plot_nn_dist('d++', d_eq, 'microcebus_nn_dist_deq')
plot_nn_dist('d+-', d_op, 'microcebus_nn_dist_dop')
plot_count_sd(all_areas, all_n,'microcebus_count_sd')

load macaque_Angelucci.mat

plot_nn_dist('d', d, 'macaque_nn_dist_d')
plot_nn_dist('d++', d_eq, 'macaque_nn_dist_deq')
plot_nn_dist('d+-', d_op, 'macaque_nn_dist_dop')
plot_count_sd(all_areas, all_n,'macaque_count_sd')

tmp = load('macaque_GrinvaldOkamoto.mat');
d = [d; tmp.d];
d_eq = [d_eq; tmp.d_eq];
d_op = [d_op; tmp.d_op];
all_areas = [all_areas; tmp.all_areas];
all_n = [all_n; tmp.all_n];

plot_nn_dist('d', d, 'macaqueAll_nn_dist_d')
plot_nn_dist('d++', d_eq, 'macaqueAll_nn_dist_deq')
plot_nn_dist('d+-', d_op, 'macaqueAll_nn_dist_dop')
plot_count_sd(all_areas, all_n,'macaqueAll_count_sd')