%% Prep
clc, clear, close all
setDefaultFigureProperties()

%% Curve Generation
% Parameter
samples = 300;
t = linspace(0, 2*pi, samples);  % parametric variable

% Parameters to control shape
a = 10;      % major amplitude (horizontal)
b = 5;    % vertical amplitude
c = 3;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = a * sin(t);
y = b * sin(2*t);
z = c * cos(t);

% Plot
figure(1);

t_s = sym("t_s","real");
xs = a * sin(t_s);
ys = b * sin(2 * t_s);
zs = c * cos(t_s);

r = [xs, ys, zs]';
r_d = simplify(diff(r, t_s));

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
nn = double(subs(ns, t_s, start));
bb = double(subs(bs, t_s, start));

plot3(x, y, z, 'k', 'LineWidth', 2);
axis equal;
xlabel('East'); ylabel('North'); zlabel('Up');
title('3D Infinity');
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
    
    kappa_arr = [kappa_arr, dtds * norm(double(subs(ts_d, t_s, ts_e)))];
    tau_arr = [tau_arr, -dtds * dot(double(subs(bs_d, t_s, ts_e)), ...
        double(subs(ns, t_s, ts_e)))];

    s_arr = [s_arr, s_arr(end) + norm(p - p_prev)];
    
    % Substitude parameter to obtain TNB vectors symbolically
    try
        t = double(subs(ts, t_s, ts_e));
        n = double(subs(ns, t_s, ts_e));
        b = double(subs(bs, t_s, ts_e));
    catch
        err_idx = [err_idx, ts_e];
        plot3(p(1), p(2), p(3), 'o', 'MarkerSize', 5);
    end

    % Create TNB Frame from local tangent and derivatives of it
    pt = [p, p + 0.5 * t];
    pn = [p, p + 0.5 * n];
    pb = [p, p + 0.5 * b];

    tt_arr = [tt_arr, t];
    nn_arr = [nn_arr, n];
    bb_arr = [bb_arr, b];

    [yaw, pitch, roll] = dcm2angle( [t,n,b]', 'zyx', 'robust');

    yaw_arr = [yaw_arr, yaw];
    pitch_arr = [pitch_arr, pitch];
    roll_arr = [roll_arr, roll];

    % Frenet-Serret
    tt = tt + ds * nn * kappa_arr(end);
    nn = nn + ds * (-tt * kappa_arr(end) + ...
        bb * tau_arr(end));
    bb = bb - ds * nn * tau_arr(end);
    
    % Create TNB Frame with Frenet-Serret
    ptt = [p, p + tt];
    pnn = [p, p + nn];
    pbb = [p, p + bb];

    if i == 1
        pt = [p, p + 1.5 * t];
        pn = [p, p + 1.5 * n];
        pb = [p, p + 1.5 * b];
        plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 4, 'Color', 'r');
        plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 4, 'Color', 'g');
        plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 4, 'Color', 'b');
        [spx, spy, spz] = sphere(10);
        surf(p(1) + 0.5*spx, p(2) + 0.5*spy, p(3) + 0.5*spz, 'FaceAlpha',0.5, 'FaceColor', 'r');
    elseif mod(i,3) == 0
        plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 2, 'Color', 'r');
        plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 2, 'Color', 'g');
        plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 2, 'Color', 'b');
    end
    % plot3(ptt(1,:), ptt(2,:), ptt(3,:), 'LineWidth', 2, 'Color', 'r', ...
    %     'LineStyle', '-.');
    % plot3(pnn(1,:), pnn(2,:), pnn(3,:), 'LineWidth', 2, 'Color', 'g', ...
    %     'LineStyle', '-.');
    % plot3(pbb(1,:), pbb(2,:), pbb(3,:), 'LineWidth', 2, 'Color', 'b', ...
    %     'LineStyle', '-.');
    pause(0.005);
    p_prev = p;
end
s_arr = s_arr(2:end);

figure(2)
plot(s_arr, kappa_arr, 'LineWidth', 4);
hold all
plot(s_arr, tau_arr, 'LineWidth', 4);
xlabel('Curvilinear Distance [m]')
ylabel('Value [1/m]')
legend('Curvature','Torsion')
title('Arclength -vs- Curvature and Torsion Profile of The Trefoil')
grid minor

display(strcat('Problematic points are', 20, num2str(err_idx')))

save infinity3d s_arr x_arr y_arr z_arr kappa_arr tau_arr tt_arr ...
    nn_arr bb_arr yaw_arr pitch_arr roll_arr