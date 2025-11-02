%% Prep
% clc, clear, close all
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
x = a * (sin(t) + 2 * sin(2 * t));
y = b * (cos(t) - 2 * cos(2 * t));
z = c * ( - sin(3 * t));

% Plot
figure(2);

t_s = sym("t_s","real");
xs = a * (sin(t_s) + 2 * sin(2 * t_s));
ys = b * (cos(t_s) - 2 * cos(2 * t_s));
zs = c * ( - sin(3 * t_s));

r = [xs, ys, zs]';
r_d = simplify(diff(r, t_s));
r_dd = simplify(diff(r_d, t_s));

ts = simplify(r_d / norm(r_d));

start = eps;
err_idx = [];
x_arr = [];
y_arr = [];
z_arr = [];
kap1_arr = [];
kap2_arr = [];
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
points_arr = [];
for i = 1:1:samples
    % cla;
    % plot3(x, y, z, 'k', 'LineWidth', 2);

    ts_e = i * 2 * pi / samples + start;
    dt = 2 * pi / samples;
    p = double(subs(r, t_s, ts_e));
    ds = norm(p - p_prev);
    dtds =  dt/ds;

    % Substitutions for the points, curvature, torsion, arclength
    x_arr = [x_arr, p(1)];
    y_arr = [y_arr, p(2)];
    z_arr = [z_arr, p(3)];
    s_arr = [s_arr, s_arr(end) + norm(p - p_prev)];

    % No torsion angular velocity
    omega = cross(double(subs(r_d, t_s, ts_e)), ...
        double(subs(r_dd, t_s, ts_e))) / norm(double(subs(r_d, t_s, ts_e)))^2;
    % dcm_prime = skew(omega) * [tt, nn, bb];
    
    omega_p = [tt, nn, bb]' * omega;
    omega_p(1) = 0;

    kap1_arr = [kap1_arr, dtds * omega_p(2,1)];
    kap2_arr = [kap2_arr, dtds * omega_p(3,1)];

    % dcm_prime = [tt, nn, bb] * skew([tt, nn, bb]' * omega);
    dcm_prime = skew(omega) * [tt, nn, bb];
    tt = tt + dt * dcm_prime(:,1);
    nn = nn + dt * dcm_prime(:,2);
    bb = bb + dt * dcm_prime(:,3);

    [yaw, pitch, roll] = dcm2angle( [tt,nn,bb]', 'zyx', 'robust');
    yaw_arr = [yaw_arr, yaw];
    pitch_arr = [pitch_arr, pitch];
    roll_arr = [roll_arr, roll];
    tt_arr = [tt_arr, tt/norm(tt)];
    nn_arr = [nn_arr, nn/norm(nn)];
    bb_arr = [bb_arr, bb/norm(bb)];

    % Create TNB Frame with Frenet-Serret
    pt = [p, p + tt];
    pn = [p, p + nn];
    pb = [p, p + bb];
    omegab = [p, p + omega];

    if mod(i,1) == 0
        % plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 2, 'Color', 'r');
        % plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 2, 'Color', 'g');
        % plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 2, 'Color', 'b');
        if i == 1
            points = plotCircle3DNew(p', ...
            [tt'],0.65,'g', 4, 1);
        elseif i == samples
            points = plotCircle3DNew(p', ...
            [tt'],0.65,'r', 4, 1);
        elseif mod(i,1) == 0
            points = plotCircle3DNew(p', ...
            [tt'],0.65,'k', 0.2, 1);
        end
        points_arr(:,:,i) = points;
        % plot3(omegab(1,:), omegab(2,:), omegab(3,:), 'LineWidth', 2, 'Color', 'm');
    end
    pause(0.005);
    p_prev = p;
end
s_arr = s_arr(2:end);

s = surf(squeeze(points_arr(1,:,:)), ...
    squeeze(points_arr(2,:,:)), ...
    squeeze(points_arr(3,:,:)));
s.FaceAlpha = 1;
s.EdgeColor = 'none';
s.FaceColor = [0.5, 0.5, 0.5];

% p_prev = double(subs(r, t_s, start));
% for i = 1:1:samples
%     ts_e = i * 2 * pi / samples + start;
%     dt = 2 * pi / samples;
%     p = double(subs(r, t_s, ts_e));
%     ds = norm(p - p_prev);
%     dtds =  dt/ds;
%     % No torsion angular velocity
%     omega = cross(double(subs(r_d, t_s, ts_e)), ...
%         double(subs(r_dd, t_s, ts_e))) / norm(double(subs(r_d, t_s, ts_e)))^2;
%     % dcm_prime = skew(omega) * [tt, nn, bb];
%     tt = tt_arr(:,i);
%     nn = nn_arr(:,i);
%     bb = bb_arr(:,i);
%     omega_p = [tt, nn, bb]' * omega;
%     kap1_arr = [kap1_arr, dtds * omega_p(2,1)];
%     kap2_arr = [kap2_arr, dtds * omega_p(3,1)];
%     p_prev = p;
% end

figure(1)
plot(s_arr, kap1_arr, 'LineWidth', 4);
hold all
plot(s_arr, kap2_arr, 'LineWidth', 4);
xlabel('Curvilinear Distance [m]')
ylabel('Value [1/m]')
legend('$$\kappa_1$$','$$\kappa_2$$','interpreter','latex')
title('Arclength -vs- Generalized Curvature Profiles of The Curve')
grid minor

display(strcat('Problematic points are', 20, num2str(err_idx')))

save trefoil_ptf s_arr x_arr y_arr z_arr kap1_arr kap2_arr tt_arr ...
    nn_arr bb_arr yaw_arr pitch_arr roll_arr