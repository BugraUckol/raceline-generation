%% Prep
clear
setDefaultFigureProperties()

%% Curve
% System
mass = 0.250;
J_xx = 10e-2;
J_yy = J_xx;
J_zz = 2*J_xx;
J_b = [J_xx, 0, 0; 0, J_yy, 0; 0, 0, J_zz];

% Curve Parameter
samples = 5000;
t_scale = 1;
t_last = t_scale*2*pi;
t_arr = linspace(0, t_last, samples);  % parametric variable
dt = t_last / (samples - 1);

% Parameters to control shape
a = 10;      % major amplitude (horizontal)
b = 4;    % vertical amplitude
c = 0.4;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = a * sin(t_arr/t_scale) + b * sin(2*t_arr/t_scale) / 5;
y = b * cos(t_arr/t_scale) - a * cos(4*t_arr/t_scale) / 8;
z = c * cos(8 * t_arr/t_scale);

g = [0,0,-9.81]';

t_s = sym("t_s","real");
xs = a * sin(t_s / t_scale) + b * sin(2 * t_s / t_scale) / 5;
ys = b * cos(t_s / t_scale) - a * cos(4 * t_s / t_scale) / 8;
zs = c * cos(8 * t_s / t_scale);

%% Calculating Required Derivatives
r = [xs, ys, zs]';
r_d = simplify(diff(r, t_s));
r_dd = simplify(diff(r_d, t_s));
r_ddd = simplify(diff(r_dd, t_s));

F_i = mass * (r_dd - g);
F_i_norm = simplify(norm(F_i));
F_i_unit = simplify(F_i / F_i_norm);

F_i_d = simplify(diff(F_i, t_s));
F_i_norm_d = simplify(diff(F_i_norm, t_s));
omega_bi_i = simplify(cross(F_i_unit,(F_i_d * F_i_norm - F_i_norm_d * F_i) / (F_i_norm^2)));
alpha_bi_i = simplify(diff(omega_bi_i, t_s));

% Initial PTF
bb = double(subs(F_i_unit, t_s, t_arr(1)));
null_tt = null(bb'/norm(bb));
nn = null_tt(:,2);
tt = null_tt(:,1);
if abs(max(cross(tt,nn) - bb)) > 0.1
    nn = null_tt(:,1);
    tt = null_tt(:,2);
end

tt_arr = [tt/norm(tt)];
nn_arr = [nn/norm(nn)];
bb_arr = [bb/norm(bb)];
c_res = [tt, nn, bb];

[a,b,c] = dcm2angle([tt, nn, bb]');
E_res = [c;b;a];

% Array Substitutions
omega_bi_i_arr = double(subs(omega_bi_i, t_s, t_arr));
alpha_bi_i_arr = double(subs(alpha_bi_i, t_s, t_arr));

alpha_bi_b = double([tt, nn, bb]' * subs(alpha_bi_i, t_s, t_arr(1)));

r_d_arr = double(subs(r_d, t_s, t_arr));
F_i_arr = double(subs(F_i, t_s, t_arr));

figure(3)
subplot(6,1,1)
plot(t_arr, F_i_arr(1,:)); hold on
subplot(6,1,2)
plot(t_arr, F_i_arr(2,:)); hold on
subplot(6,1,3)
plot(t_arr, F_i_arr(3,:)); hold on
subplot(6,1,4)
plot(t_arr, alpha_bi_i_arr(1,:)); hold on
subplot(6,1,5)
plot(t_arr, alpha_bi_i_arr(2,:)); hold on
subplot(6,1,6)
plot(t_arr, alpha_bi_i_arr(3,:)); hold on
figure(4)
subplot(3,1,1)
plot(t_arr, omega_bi_i_arr(1,:)); hold on
subplot(3,1,2)
plot(t_arr, omega_bi_i_arr(2,:)); hold on
subplot(3,1,3)
plot(t_arr, omega_bi_i_arr(3,:)); hold on

%% Extracting Differentially Flat Inputs
% Initializing Differential Flatness Arrays
df_omega_bi_b_arr = [tt, nn, bb]' * omega_bi_i_arr(:,1);
df_vel_bi_b_arr = [tt, nn, bb]' * r_d_arr(:,1);
df_eul_arr = E_res;
df_thrust_arr = [[tt, nn, bb]' * F_i_arr(:,1)];
df_moment_arr = [J_b * alpha_bi_b + cross(df_omega_bi_b_arr(:,1), J_b * df_omega_bi_b_arr(:,1))];

% Loop
counter = 1;
for t = t_arr(2:end)
    counter = counter + 1;
    omega_prev = c_res' * omega_bi_i_arr(:,counter-1);

    if norm(omega_prev) > eps
        nk = omega_prev / norm(omega_prev);
        dk = norm(omega_prev) * dt;
        R = eye(3) + skew(nk) * sin(dk) + skew(nk) * skew(nk) * (1 - cos(dk));
    else
        R = eye(3);
    end

    c_res = c_res * R;

    tt = c_res(:,1);
    nn = c_res(:,2);
    bb = c_res(:,3);

    tt_arr = [tt_arr, tt];
    nn_arr = [nn_arr, nn];
    bb_arr = [bb_arr, bb];

    [a,b,c] = dcm2angle([tt, nn, bb]');
    E_res = [c;b;a];
    df_eul_arr = [df_eul_arr, E_res];

    omega_bi_b = [tt, nn, bb]' * omega_bi_i_arr(:,counter);
    omega_bi_b(3) = 0;
    df_omega_bi_b_arr = [df_omega_bi_b_arr, omega_bi_b];
    df_vel_bi_b_arr = [df_vel_bi_b_arr, [tt, nn, bb]' * r_d_arr(:,counter)];
    alpha_bi_b = [tt, nn, bb]' * alpha_bi_i_arr(:,counter);

    df_thrust_arr = [df_thrust_arr, [0;0;[0,0,1] * [tt, nn, bb]' * F_i_arr(:,counter)]];
    df_moment_arr = [df_moment_arr, (J_b * alpha_bi_b + cross(omega_bi_b, J_b * omega_bi_b))];
end

figure(2)
subplot(4,1,1)
plot(t_arr, df_thrust_arr(3,:)); hold on
subplot(4,1,2)
plot(t_arr, df_moment_arr(1,:)); hold on
subplot(4,1,3)
plot(t_arr, df_moment_arr(2,:)); hold on
subplot(4,1,4)
plot(t_arr, df_moment_arr(3,:)); hold on

%% Differential Flatness Test by Open Loop Simulation
tt = tt_arr(:,1);
nn = nn_arr(:,1);
bb = bb_arr(:,1);
P_res = [x(1), y(1), z(1)]';
P_res_all = P_res;
v_res = df_vel_bi_b_arr(:,1);
w_res = df_omega_bi_b_arr(:,1);
c_res = [tt, nn, bb];
[a,b,c] = dcm2angle([tt, nn, bb]');
E_res = [c;b;a];

% close all
figure(1)
subplot(1,2,1)
plot3(x,y,z); hold on; daspect([1,1,1])
bx = [P_res, P_res + c_res(:,1)];
by = [P_res, P_res + c_res(:,2)];
bz = [P_res, P_res + c_res(:,3)];
plot3(bx(1,:), bx(2,:), bx(3,:), 'LineWidth', 2, 'Color', 'r');
plot3(by(1,:), by(2,:), by(3,:), 'LineWidth', 2, 'Color', 'g');
plot3(bz(1,:), bz(2,:), bz(3,:), 'LineWidth', 2, 'Color', 'b');
for k=2:1:samples
    T_res = df_thrust_arr(:,k-1);
    F_res = T_res + mass * c_res' * g;
    M_res = df_moment_arr(:,k-1);

    P_res = P_res + dt * c_res * v_res;
    v_res = v_res + dt * (F_res/mass - cross(w_res, v_res));
    w_current = w_res;
    if norm(w_current) > eps
        nk = w_current / norm(w_current);
        dk = norm(w_current) * dt;
        R = eye(3) + skew(nk) * sin(dk) + skew(nk) * skew(nk) * (1 - cos(dk));
    else
        R = eye(3);
    end

    c_res = c_res * R;
    w_res = w_res + dt * J_b^-1 * (M_res - cross(w_res, J_b * w_res));
    P_res_all = [P_res_all, P_res];

    bx = [P_res, P_res + c_res(:,1)];
    by = [P_res, P_res + c_res(:,2)];
    bz = [P_res, P_res + c_res(:,3)];

    v_res
    df_vel_bi_b_arr(:,k)
    w_res
    df_omega_bi_b_arr(:,k)
    
    % if mod(k,400) == 0
    % % plot3(bx(1,:), bx(2,:), bx(3,:), 'LineWidth', 2, 'Color', 'r');
    % % plot3(by(1,:), by(2,:), by(3,:), 'LineWidth', 2, 'Color', 'g');
    % % plot3(bz(1,:), bz(2,:), bz(3,:), 'LineWidth', 2, 'Color', 'b');
    % associated_r = subs(r,t_s,t_arr(k));
    % inter_vecs = double([P_res, associated_r]);
    % plot3(inter_vecs(1,:), inter_vecs(2,:), inter_vecs(3,:),'k','LineWidth',0.7)
    % end
end
plot3(P_res_all(1,1:10:end), P_res_all(2,1:10:end), P_res_all(3,1:10:end),...
        'LineWidth',2,'Color', 'm')

%% Differential Flatness Test by Open Loop Kinematic Simulation
subplot(1,2,2)
tt = tt_arr(:,1);
nn = nn_arr(:,1);
bb = bb_arr(:,1);
c_res = [tt, nn, bb];
P_res = [x(1), y(1), z(1)]';
P_res_all = P_res;
v_res_arr = double(subs(r_d, t_s, t_arr));
w_res_arr = double(subs(omega_bi_i, t_s, t_arr));
w_res(3) = 0;
figure(1)
plot3(x,y,z); hold on; daspect([1,1,1])
bx = [P_res, P_res + c_res(:,1)];
by = [P_res, P_res + c_res(:,2)];
bz = [P_res, P_res + c_res(:,3)];

plot3(bx(1,:), bx(2,:), bx(3,:), 'LineWidth', 2, 'Color', 'r');
plot3(by(1,:), by(2,:), by(3,:), 'LineWidth', 2, 'Color', 'g');
plot3(bz(1,:), bz(2,:), bz(3,:), 'LineWidth', 2, 'Color', 'b');
for k=1:1:samples
    P_res = P_res + dt * v_res_arr(:,k);
    v_res = c_res' * v_res_arr(:,k);

    w_current = c_res' * w_res_arr(:,k);
    if norm(w_current) > eps
        nk = w_current / norm(w_current);
        dk = norm(w_current) * dt;
        R = eye(3) + skew(nk) * sin(dk) + skew(nk) * skew(nk) * (1 - cos(dk));
    else
        R = eye(3);
    end

    c_res = c_res * R;
    w_res = c_res' * w_res_arr(:,k);
    w_res(3) = 0;
    w_res;

    if mod(k,400) == 0

        bx = [P_res, P_res + c_res(:,1) * 0.3];
        by = [P_res, P_res + c_res(:,2) * 0.3];
        bz = [P_res, P_res + c_res(:,3) * 0.3];
        plot3(bx(1,:), bx(2,:), bx(3,:), 'LineWidth', 2, 'Color', 'r');
        plot3(by(1,:), by(2,:), by(3,:), 'LineWidth', 2, 'Color', 'g');
        plot3(bz(1,:), bz(2,:), bz(3,:), 'LineWidth', 2, 'Color', 'b');
        % associated_r = subs(r,t_s,t_arr(k));
        % inter_vecs = double([P_res, associated_r]);
        % plot3(inter_vecs(1,:), inter_vecs(2,:), inter_vecs(3,:),'k','LineWidth',0.7)


    end
    P_res_all = [P_res_all, P_res];
    plot3(P_res_all(1,end-1:end), P_res_all(2,end-1:end), P_res_all(3,end-1:end),...
        'LineWidth',2,'Color', 'm')
end

%% Figure2
r_dd_arr = double(subs(r_dd, t_s, t_arr));
tot_acc_arr = r_dd_arr - g;
figure(5)
subplot(2,2,1)
plot3(x,y,z);
subplot(2,2,2)
plot3(r_d_arr(1,:),r_d_arr(2,:),r_d_arr(3,:)); 
subplot(2,2,3)
plot3(r_dd_arr(1,:),r_dd_arr(2,:),r_dd_arr(3,:)); 
subplot(2,2,4)
plot3(tot_acc_arr(1,:),tot_acc_arr(2,:),tot_acc_arr(3,:)); 

%% Save
% df_thrust_arr = df_thrust_arr(3,:);
% df_moment_arr = df_moment_arr;
% save diff_flat_generic_curve samples dt t_arr df_thrust_arr df_moment_arr ...
%     df_vel_bi_b_arr df_omega_bi_b_arr df_eul_arr

















