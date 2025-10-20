%% Prep
clear
setDefaultFigureProperties()

%% Curve and Differential Flatness Generation
% System
mass = 0.250;
J_xx = 10e-2;
J_yy = J_xx;
J_zz = J_xx;
J_b = [J_xx, 0, 0; 0, J_yy, 0; 0, 0, J_zz];

% Curve Parameter
samples = 30000;
t_last = 2*pi;
t_arr = linspace(0, t_last, samples);  % parametric variable
dt = t_last / (samples - 1);

% Parameters to control shape
a = 10;      % major amplitude (horizontal)
b = 4;    % vertical amplitude
c = 0.4;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = a * sin(t_arr) + b * sin(2*t_arr) / 5;
y = b * cos(t_arr) - a * cos(4*t_arr) / 8;
z = c * cos(8 * t_arr);

g = [0,0,-9.81]';

% Plot
figure(10);

t_s = sym("t_s","real");
xs = a * sin(t_s) + b * sin(2 * t_s) / 5;
ys = b * cos(t_s) - a * cos(4 * t_s) / 8;
zs = c * cos(8 * t_s);

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
% omega_bi_i2 = simplify(cross(F_i_unit, simplify(diff(F_i_unit, t_s))))

simplify(dot(omega_bi_i,F_i_unit)); % Must be zero since no angular velocity is required on this axis

alpha_bi_i = simplify(diff(omega_bi_i, t_s));

tt = double(subs(F_i_unit, t_s, t_arr(1)));
null_tt = null(tt');
nn = null_tt(:,1);
bb = null_tt(:,2);
cross(tt,nn) - bb

tt_arr = [tt];
nn_arr = [nn];
bb_arr = [bb];

omega_bi_b = [tt, nn, bb]' * subs(omega_bi_i, t_s, t_arr(1));
alpha_bi_b = [tt, nn, bb]' * subs(alpha_bi_i, t_s, t_arr(1));

Thrust_arr = [double([tt, nn, bb]' * subs(F_i, t_s, t_arr(1)))];
Moment_arr = [double(J_b * alpha_bi_b + cross(omega_bi_b, J_b * omega_bi_b))];
Uvw_arr = [double([tt, nn, bb]' * subs(r_d, t_s, t_arr(1)))];
Omega_arr = [double([tt, nn, bb]' * subs(omega_bi_i, t_s, t_arr(1)))];

c_res = [tt, nn, bb];
for t = t_arr(2:end)
    omega = c_res' * double(subs(omega_bi_i, t_s, t));
    
    if norm(omega) > eps
    nk = omega / norm(omega);
    dk = norm(omega) * dt;
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

    omega_bi_b = [tt, nn, bb]' * double(subs(omega_bi_i, t_s, t));
    omega_bi_b(1) = 0;
    alpha_bi_b = [tt, nn, bb]' * double(subs(alpha_bi_i, t_s, t));

    Thrust_arr = [Thrust_arr, double([tt, nn, bb]' * subs(F_i, t_s, t))];
    Moment_arr = [Moment_arr, double(J_b * alpha_bi_b + cross(omega_bi_b, J_b * omega_bi_b))];
    Uvw_arr = [Uvw_arr, double([tt, nn, bb]' * subs(r_d, t_s, t))];
    Omega_arr = [Omega_arr, double([tt, nn, bb]' * subs(omega_bi_i, t_s, t))];
end

%% Differential Flatness Test by Open Loop Simulation
close all
figure(1)
subplot(1,2,1)
tt = tt_arr(:,1);
nn = nn_arr(:,1);
bb = bb_arr(:,1);
P_res = [x(1), y(1), z(1)]';
P_res_all = P_res;
v_res = [tt, nn, bb]' * double(subs(r_d, t_s, t_arr(1)));
w_res = [tt, nn, bb]' * double(subs(omega_bi_i, t_s, t_arr(1)));
w_res(1) = 0;
c_res = [tt, nn, bb];
figure(1)
plot3(x,y,z); hold on; daspect([1,1,1])
bx = [P_res, P_res + c_res(:,1)];
by = [P_res, P_res + c_res(:,2)];
bz = [P_res, P_res + c_res(:,3)];

plot3(bx(1,:), bx(2,:), bx(3,:), 'LineWidth', 2, 'Color', 'r');
plot3(by(1,:), by(2,:), by(3,:), 'LineWidth', 2, 'Color', 'g');
plot3(bz(1,:), bz(2,:), bz(3,:), 'LineWidth', 2, 'Color', 'b');

for k=1:1:samples
    T_res = Thrust_arr(:,k);
    F_res = T_res;
    M_res = Moment_arr(:,k);

    P_res = P_res + dt * c_res * v_res;
    v_res = v_res + dt * (F_res/mass + c_res'*[0;0;g(3)] - cross(w_res, v_res));
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
    w_res(1) = 0;
    w_res;
    P_res_all = [P_res_all, P_res];

    bx = [P_res, P_res + c_res(:,1)];
    by = [P_res, P_res + c_res(:,2)];
    bz = [P_res, P_res + c_res(:,3)];
    
    if mod(k,400) == 0
    % plot3(bx(1,:), bx(2,:), bx(3,:), 'LineWidth', 2, 'Color', 'r');
    % plot3(by(1,:), by(2,:), by(3,:), 'LineWidth', 2, 'Color', 'g');
    % plot3(bz(1,:), bz(2,:), bz(3,:), 'LineWidth', 2, 'Color', 'b');
    associated_r = subs(r,t_s,t_arr(k));
    inter_vecs = double([P_res, associated_r]);
    plot3(inter_vecs(1,:), inter_vecs(2,:), inter_vecs(3,:),'k','LineWidth',0.7)
    end
    plot3(P_res_all(1,end-1:end), P_res_all(2,end-1:end), P_res_all(3,end-1:end),...
        'LineWidth',2,'Color', 'm')
end

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
w_res(1) = 0;
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
    w_res(1) = 0;
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

%% Save
save generic_path_30_000_steps_diff_flat dt Thrust_arr Moment_arr ...
    Uvw_arr Omega_arr

















