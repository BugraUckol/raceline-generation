clc, clear, close all

% Parameter
t = linspace(0, 2*pi, 1000);  % parametric variable

% Parameters to control shape
a = 10;      % major amplitude (horizontal)
b = 5;    % vertical amplitude
c = 3;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = a * sin(t);
y = b * sin(2*t);
z = c * cos(t);

% Plot
figure;
plot3(x, y, z, 'k', 'LineWidth', 2);
    grid on; axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
    title('3D Infinity Curve (Non-Intersecting)');
    view(135, 30); % adjust view angle for clarity
    hold on

t_s = sym("t_s","real");

xs = a * sin(t_s);
ys = b * sin(2*t_s);
zs = c * cos(t_s);

r = [xs, ys, zs]';
r_d = diff(r, t_s);

ts = r_d / norm(r_d);

ts_d = simplify(diff(ts, t_s));

ns = simplify(ts_d / norm(ts_d));

bs = simplify(cross(ts, ns));

for i = -10000:5:10000
    cla;
    plot3(x, y, z, 'k', 'LineWidth', 2);

    ts_e = i * 2*pi / 1000;
    p = double(subs(r, t_s, ts_e));
    
    try
        tic
        t = double(subs(ts, t_s, ts_e));
        n = double(subs(ns, t_s, ts_e));
        b = double(subs(bs, t_s, ts_e));
        toc
    end
    
    pt = [p, p + t];
    pn = [p, p + n];
    pb = [p, p + b];

    plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 2, 'Color', 'r');
    plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 2, 'Color', 'g');
    plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 2, 'Color', 'b');
    drawnow
    toc
end
