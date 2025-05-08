%% Prep
clc, clear, close all

%% Curve Generation
% Parameter
t = linspace(0, 2*pi, 100);  % parametric variable

% Parameters to control shape
a = 10;      % major amplitude (horizontal)
b = 5;    % vertical amplitude
c = 3;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = 0.1 * t .* a .* sin(t);
y = 0.1 * t .* b .* sin(2*t);
z = 0.1 * t .* c .* cos(t);

% Plot
figure(1);

syms t_s
t_s = sym("t_s","real");
xs = 0.1 * t_s * a * sin(t_s);
ys = 0.1 * t_s * b * sin(2*t_s);
zs = 0.1 * t_s * c * cos(t_s);

r = [xs, ys, zs]';
r_d = diff(r, t_s);

ts = r_d / norm(r_d);

ts_d = diff(ts, t_s);

ns = ts_d / norm(ts_d);

bs = cross(ts, ns);

bs_d = diff(bs, t_s);

start = eps;
err_idx = [];
x_arr = [];
y_arr = [];
z_arr = [];
kappa_arr = [];
tau_arr = [];
s_arr = [0];
p_prev = double(subs(r, t_s, start));
tt = double(subs(ts, t_s, start));
nn = double(subs(ns, t_s, start));
bb = double(subs(bs, t_s, start));
for i = 1:100
    clf;
    plot3(x, y, z, 'k', 'LineWidth', 2);
    grid on; axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
    title('3D Infinity Curve (Non-Intersecting)');
    view(135, 30); % adjust view angle for clarity
    hold on

    ts_e = i * 2*pi / 1000 + start;
    p = double(subs(r, t_s, ts_e));

    % Substitutions for the points, curvature, torsion, arclength
    dtds = 2*pi / 1000 / norm(p - p_prev);
    x_arr = [x_arr, p(1)];
    y_arr = [y_arr, p(2)];
    z_arr = [z_arr, p(3)];
    kappa_arr = [kappa_arr, dtds * norm(double(subs(ts_d, t_s, ts_e)))];
    tau_arr = [tau_arr, -dtds * norm(double(subs(bs_d, t_s, ts_e)))];
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

    % Frenet-Serret
    tt = tt + norm(p - p_prev) * nn * kappa_arr(end);
    nn = nn + norm(p - p_prev) * (-tt * kappa_arr(end) + ...
        bb * tau_arr(end));
    bb = bb - norm(p - p_prev) * nn * tau_arr(end);
    
    % Create TNB Frame with Frenet-Serret
    ptt = [p, p + tt];
    pnn = [p, p + nn];
    pbb = [p, p + bb];

    plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 2, 'Color', 'r');
    plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 2, 'Color', 'g');
    plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 2, 'Color', 'b');

    plot3(ptt(1,:), ptt(2,:), ptt(3,:), 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-.');
    plot3(pnn(1,:), pnn(2,:), pnn(3,:), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '-.');
    plot3(pbb(1,:), pbb(2,:), pbb(3,:), 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-.');
    pause(0.005);
    p_prev = p;
end
s_arr = s_arr(2:end);

figure(2)
plot(s_arr, kappa_arr, 'LineWidth', 2);
hold all
plot(s_arr, tau_arr, 'LineWidth', 2);

display(strcat('Problematic points are', 20, num2str(err_idx')))

save three_d_infinity x y z