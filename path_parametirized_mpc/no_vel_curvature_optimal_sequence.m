%% Prep
clc, clear%, close all
% Clear is especially important since the opti object should not be 
%   given constraints from the past

%% Casadi Imports
addpath("/Users/bugrauckol/Documents/share/casadi-3")
import casadi.*

%% Initial conditions and constants
t0 = 0;
v0 = 1.0;
e_y0 = 0.0;
e_p0 = 0.0;

lw = 2.971; % Wheelbase Indy Autonomous
steering_lim = 0.5;
max_lat_acc = 1*9.81;

%% Load
load yas_marina_8000.mat
% load yas_marina_8000_reparam.mat
x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr = spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr;
% d = designfilt("lowpassfir", ...
%     PassbandFrequency=0.05,StopbandFrequency=0.2, ...
%     PassbandRipple=1,StopbandAttenuation=60, ...
%     DesignMethod="equiripple");
% x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:, 4) = ...
%     filtfilt(d, x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:, 4));
% figure
% hold on
% plot(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,4))
% plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,4))

%% Setting Optimization Problem
size_vec = size(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr);
lim = size_vec(1) - 1;
N = lim; %size_vec(1) - 1;
opti = casadi.Opti(); % Optimization problem

% ---- decision variables ---------
X = opti.variable(3, N+1); % state trajectory in path frame
t = X(1,:);
ey = X(2,:);
ep = X(3,:);

U = opti.variable(1,N);   %steering

% ---- objective          ---------
% opti.minimize(t(end)); % minimize time
opti.minimize(1.0 * U(1,:) * U(1,:)' + 0.0 * t(end)); % minimize steering

% ---- dynamic constraints --------
% x' = [t, ey, ep]
f = @(x,u,kappa) [
    (1 - kappa * x(2)) / (v0 * cos(x(3)));
    (v0 * sin(x(3))) * ((1 - kappa * x(2)) / (v0 * cos(x(3))));
    (v0 * tan(u(1)) / lw) * ((1 - kappa * x(2)) / (v0 * cos(x(3)))) - kappa; 
   ];

% guess_steering = [];
% guess_speed = [];

for k=1:N % loop over control intervals
   kappa = x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k,4);
   ds = x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k + 1, 5) ...
       - x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k, 5);

   % guess_steering = [guess_steering, tan(kappa * lw)];
   % guess_speed = [guess_speed, sqrt(max_lat_acc / abs(tan(guess_steering(end)) / lw))];
   
   % % Runge-Kutta 4 integration
   % k1 = f(X(:,k),         U(:,k), kappa);
   % k2 = f(X(:,k)+ds/2*k1, U(:,k), kappa);
   % k3 = f(X(:,k)+ds/2*k2, U(:,k), kappa);
   % k4 = f(X(:,k)+ds*k3,   U(:,k), kappa);
   % x_next = X(:,k) + ds/6*(k1+2*k2+2*k3+k4);
   
   % 1st Order Explicit Euler's Integration
   k1 = f(X(:,k), U(:,k), kappa);
   x_next = X(:,k) + ds * k1;
   
   % Multiple shooting
   opti.subject_to(X(:,k+1)==x_next); % close the gaps

   opti.subject_to(X(2,k+1) <= x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k + 1, 10));
   opti.subject_to(X(2,k+1) >= x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k + 1, 11));

   opti.subject_to(X(2,k+1) * x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k + 1, 4) < 0.9);

   opti.subject_to(X(1,k+1) > X(1,k));
end

% Steering is limited
opti.subject_to(-steering_lim <= U(1,:))
opti.subject_to(U(1,:) <= steering_lim);

% Angle Limits
opti.subject_to(X(3,:) <= pi/2.1);
opti.subject_to(X(3,:) >= -pi/2.1);

warn_res = warndlg('COMMENTED-OUT STEER RATE LIMIT BUT MIGHT BE NECESSARY IN THE FUTURE');
% opti.subject_to((U(1,2:end) - U(1,1:end-1)).^2 <= (0.001*steering_lim)^2);

% ---- Initial guess
% Problem is symmetric when the pendulum is pointing down or up
%   adding a random initial guess helps it to solve
% opti.set_initial(U(1,:), guess_steering);

% ---- initial conditions --------
opti.subject_to(t(1) == t0);
opti.subject_to(ey(end) == ey(1));
opti.subject_to(ep(end) == ep(1));
% opti.subject_to(ey(1) == e_y0);
% opti.subject_to(ep(1) == e_p0);

% ---- final conditions -------
% opti.subject_to(ey(end) == e_y0);
% opti.subject_to(ep(end) == e_p0);

% ---- final conditions -------
% opti.subject_to(psi(end) == 0.0);
% opti.subject_to(pos_y(end) == 0.0);
% opti.subject_to(theta(end) == 0.0);
% opti.subject_to(theta_d(end) == 0.0);

% ---- set initial guess
% opti.set_initial(U, 0.5);

%% Solve NLP            ------
opts = struct();
opts.ipopt = struct();
opts.ipopt.tol = 1e-3;  % Increase tolerance
opts.ipopt.acceptable_tol = 1e-4;  % Reduce acceptable tolerance
opts.ipopt.max_iter = 10000;
opts.ipopt.mu_init = 1e-3;  % Initial value for the barrier parameter
opti.solver('ipopt', opts); % Set numerical backend
sol = opti.solve();   % actual solve

%% Extract Input Sequence and Save
distance_arr = (0:N) * ds;
time_arr = sol.value(t);
vel_arr = sol.value(t) ./ sol.value(t) * v0;
lat_err_arr = sol.value(ey);
ang_err_arr = sol.value(ep);

distance_input_arr = (0:N-1) * ds;
steering_arr = sol.value(U(1,:));
throttle_arr = sol.value(U(1,:));

steering_ts = timeseries(steering_arr, time_arr(1:end-1));
throttle_ts = timeseries(throttle_arr, time_arr(1:end-1));

%% Transform Predictions to Map
traj_x = x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,1) + lat_err_arr' ...
    .* -sin(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,3));
traj_y = x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,2) + lat_err_arr' ...
    .* cos(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,3));

%% Post Process
set(0,'DefaultFigureWindowStyle','docked')

figure(1)
subplot(2,2,1)
hold on
plot(distance_arr, lat_err_arr, 'LineWidth', 2)
plot(distance_arr, x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1, 10), '-.', 'LineWidth', 2)
plot(distance_arr, x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1, 11), '-.', 'LineWidth', 2)
title('Distance vs Lat. Error')
grid minor
subplot(2,2,2)
hold on
plot(distance_arr, ang_err_arr, 'LineWidth', 2)
title('Distance vs Ang. Error')
grid minor
subplot(2,2,3)
hold on
plot(distance_arr, vel_arr, 'LineWidth', 2)
plot(distance_arr(1:end-1), throttle_arr, 'LineWidth', 2)
title('Vel')
grid minor
subplot(2,2,4)
hold on
plot(distance_arr, time_arr, 'LineWidth', 2)
title('Distance vs. Time')
grid minor

figure(2)
hold on
plot(time_arr(1:end-1), sol.value(steering_arr), 'LineWidth', 2);
legend('$\delta$', 'Location', 'northeast', 'Interpreter', ...
    'latex')
title('Time ve Inputs')
grid minor

load yas_marina_min_curve
figure(3)
hold on
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,1), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,2), '--', 'LineWidth', ...
    0.5, 'Color', [0.05,0.05,0.05]);
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,6), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,7), 'LineWidth', 2, ...
    'Color', [0.05,0.05,0.05]);
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,8), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1:lim+1,9), 'LineWidth', 2, ...
    'Color', [0.05,0.05,0.05]);
plot(yas_marina_min_curve(:,1),yas_marina_min_curve(:,2), 'LineWidth', ...
    2, 'Color', 'cyan');
plot(traj_x(1:end-1), traj_y(1:end-1), 'LineWidth', 1.8, 'Color', ...
    'magenta');
title('Map')
legend('CenterLine','LeftMargin','RightMargin','MinCurvature', ...
    'OptimalTrajectory')
grid minor
title('Map')
daspect([1,1,1])

%% Save
save 8000_result_first_iter_long_wheelbase_rate_lim x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr ...
distance_arr time_arr vel_arr lat_err_arr ang_err_arr ...
distance_input_arr steering_arr throttle_arr steering_ts throttle_ts