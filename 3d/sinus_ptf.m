%% Prep
clc, clear, close all
setDefaultFigureProperties()

%% Curve Generation
% Parameter
samples = 1000;
t = linspace(0, 2*pi, samples);  % parametric variable

% Parameters to control shape
a = 3;      % major amplitude (horizontal)
b = 3;    % vertical amplitude
c = 5;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = t;
y = sin(t);
z = 0*t;

% Plot
figure(1);

t_s = sym("t_s","real");
xs = t_s;
ys = sin(t_s);
zs = 0;

r = [xs, ys, zs]';
r_d = simplify(diff(r, t_s));
r_dd = simplify(diff(r_d, t_s));

ts = simplify(r_d / norm(r_d));

ts_d = simplify(diff(ts, t_s));

ns = simplify(ts_d / norm(ts_d));

bs = simplify(cross(ts, ns));

bs_d = simplify(diff(bs, t_s));

start = eps;
err_idx = [];
x_arr = [];
y_arr = [];
z_arr = [];
kappa_arr = [];
tau_arr = [];
s_arr = [0];
tt_arr = [];
nn_arr = [];
bb_arr = [];
yaw_arr = [];
pitch_arr = [];
roll_arr = [];
p_prev = double(subs(r, t_s, start));
tt = double(subs(ts, t_s, start));
null_tt = null(tt');
nn = null_tt(:,1);
bb = null_tt(:,2);

plot3(x, y, z, 'k', 'LineWidth', 2);
axis equal;
xlabel('East'); ylabel('North'); zlabel('Up');
title('Trefoil');
view(135, 30); % adjust view angle for clarity
hold on
plot3(x, y, z, 'k', 'LineWidth', 2);
box on;
for i = 1:1:samples
    % cla;
    % plot3(x, y, z, 'k', 'LineWidth', 2);

    ts_e = i * 2 * pi / samples + start;
    dt = 2 * pi / samples;
    p = double(subs(r, t_s, ts_e));
    ds = norm(p - p_prev);

    % Substitutions for the points, curvature, torsion, arclength
    dtds =  dt/ds;
    x_arr = [x_arr, p(1)];
    y_arr = [y_arr, p(2)];
    z_arr = [z_arr, p(3)];

    % kappa_arr = [kappa_arr, dtds * norm(double(subs(ts_d, t_s, ts_e)))];
    % tau_arr = [tau_arr, -dtds * dot(double(subs(bs_d, t_s, ts_e)), ...
    %     double(subs(ns, t_s, ts_e)))];
    % 
    % s_arr = [s_arr, s_arr(end) + norm(p - p_prev)];
    % 
    % % Substitude parameter to obtain TNB vectors symbolically
    % try
    %     t = double(subs(ts, t_s, ts_e));
    %     n = double(subs(ns, t_s, ts_e));
    %     b = double(subs(bs, t_s, ts_e));
    % catch
    %     err_idx = [err_idx, ts_e];
    %     plot3(p(1), p(2), p(3), 'o', 'MarkerSize', 5);
    % end
    % 
    % % Create TNB Frame from local tangent and derivatives of it
    % pt = [p, p + 0.5 * t];
    % pn = [p, p + 0.5 * n];
    % pb = [p, p + 0.5 * b];
    % 
    tt_arr = [tt_arr, t];
    % nn_arr = [nn_arr, n];
    % bb_arr = [bb_arr, b];
    % 
    % [yaw, pitch, roll] = dcm2angle( [t,n,b]', 'zyx', 'robust');
    % 
    % yaw_arr = [yaw_arr, yaw];
    % pitch_arr = [pitch_arr, pitch];
    % roll_arr = [roll_arr, roll];
    % 
    % Frenet-Serret
    omega = cross(double(subs(r_d, t_s, ts_e)), ...
        double(subs(r_dd, t_s, ts_e))) / norm(double(subs(r_d, t_s, ts_e)))^2;
    % dcm_prime = skew(omega) * [tt, nn, bb];
    omega_p = [tt, nn, bb]' * omega;
    omega_p(1) = 0
    dcm_prime = [tt, nn, bb] * skew([tt, nn, bb]' * omega);
    tt = tt + dt * dcm_prime(:,1);
    nn = nn + dt * dcm_prime(:,2);
    bb = bb + dt * dcm_prime(:,3);

    % dcm_prime = skew(omega) * [tt, nn, bb];
    % tt = tt + dt * dcm_prime(:,1);
    % nn = nn + dt * dcm_prime(:,2);
    % bb = bb + dt * dcm_prime(:,3);

    % tt = t;
    % null_t = null(t');
    % nn = null_tt(:,1);
    % bb = null_tt(:,2);

    % Create TNB Frame with Frenet-Serret
    ptt = [p, p + tt];
    pnn = [p, p + nn];
    pbb = [p, p + bb];

    if i == 19999999
        pt = [p, p + 1.5 * t];
        pn = [p, p + 1.5 * n];
        pb = [p, p + 1.5 * b];
        plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 4, 'Color', 'r');
        plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 4, 'Color', 'g');
        plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 4, 'Color', 'b');
        [spx, spy, spz] = sphere(10);
        surf(p(1) + 0.5*spx, p(2) + 0.5*spy, p(3) + 0.5*spz, 'FaceAlpha',0.5, 'FaceColor', 'r');
    elseif mod(i,5) == 0
        % plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 2, 'Color', 'r');
        % plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 2, 'Color', 'g');
        % plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 2, 'Color', 'b');

        plot3(ptt(1,:), ptt(2,:), ptt(3,:), 'LineWidth', 2, 'Color', 'r', ...
            'LineStyle', '-.');
        plot3(pnn(1,:), pnn(2,:), pnn(3,:), 'LineWidth', 2, 'Color', 'g', ...
            'LineStyle', '-.');
        plot3(pbb(1,:), pbb(2,:), pbb(3,:), 'LineWidth', 2, 'Color', 'b', ...
            'LineStyle', '-.');
    end
    pause(0.005);
    p_prev = p;
end
s_arr = s_arr(2:end);

figure(2)
plot(s_arr, kappa_arr, 'LineWidth', 4);
hold on
plot(s_arr, tau_arr, 'LineWidth', 4);
xlabel('Curvilinear Distance [m]')
ylabel('Value [1/m]')
legend('Curvature','Torsion')
title('Arclength -vs- Curvature and Torsion Profile of The Trefoil')
grid minor

display(strcat('Problematic points are', 20, num2str(err_idx')))

% % % save trefoil s_arr x_arr y_arr z_arr kappa_arr tau_arr tt_arr ...
% % %     nn_arr bb_arr yaw_arr pitch_arr roll_arr