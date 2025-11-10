num_of_samples_arr = [50, 100, 200, 500, 1000, 2000];
mc_time_arr = [];
mc_distance_arr = [];

for num_of_samples = num_of_samples_arr
    alternative_ptf_generic;
    ptf_optimization;
    mc_time_arr = [mc_time_arr, opti.stats.t_wall_total];
    mc_distance_arr = [mc_distance_arr, norm(...
        [x_arr(end), y_arr(end), z_arr(end)]' - P_res)];
    clearvars -except num_of_samples num_of_samples_arr mc_time_arr ...
        mc_distance_arr
end

save simResults

figure(40)
ax = plotyy(num_of_samples_arr, mc_time_arr,num_of_samples_arr, ...
    mc_distance_arr);
xlabel('Number of Samples')
ylabel('Time [s]')
ylabel(ax(2), 'Deviation [m]')





















