diff_flats = load('generic_path_30_000_steps_diff_flat');

samples = 30000;
t_last = 2*pi;
t_arr = linspace(0, t_last, samples);

saved_samples = max(size(diff_flats.Uvw_arr));
samples = 2000;%% CMON

ratio = round(saved_samples / samples);

prevsol.N = samples;
prevsol.time_arr = t_arr(1:ratio:saved_samples);
prevsol.en_arr = zeros(samples,1);
prevsol.eb_arr = zeros(samples,1);
prevsol.ephi_arr = zeros(samples,1); %% THESE ARE NOT TO BE ZERO!!!!
prevsol.ethe_arr = zeros(samples,1); %% THESE ARE NOT TO BE ZERO!!!!
prevsol.epsi_arr = zeros(samples,1); %% THESE ARE NOT TO BE ZERO!!!!
prevsol.u_arr = diff_flats.Uvw_arr(1,1:ratio:saved_samples);
prevsol.v_arr = diff_flats.Uvw_arr(2,1:ratio:saved_samples);
prevsol.w_arr = diff_flats.Uvw_arr(3,1:ratio:saved_samples);
prevsol.p_arr = diff_flats.Omega_arr(1,1:ratio:saved_samples);
prevsol.q_arr = diff_flats.Omega_arr(2,1:ratio:saved_samples);
prevsol.r_arr = diff_flats.Omega_arr(3,1:ratio:saved_samples);
prevsol.T_com_arr = diff_flats.Thrust_arr(1,1:ratio:saved_samples);
prevsol.Mx_com_arr = diff_flats.Moment_arr(1,1:ratio:saved_samples);
prevsol.My_com_arr = diff_flats.Moment_arr(2,1:ratio:saved_samples);
prevsol.Mz_com_arr = diff_flats.Moment_arr(3,1:ratio:saved_samples);

save prevsol_diff_flat prevsol