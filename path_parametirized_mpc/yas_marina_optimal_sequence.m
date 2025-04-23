%% Prep
clc, clear, close all
% Clear is especially important since the opti object should not be 
%   given constraints from the past

%% Casadi Imports
addpath("C:\Program Files\casadi-3.6.7-windows64-matlab2018b")
import casadi.*

%% Initial conditions and constants
t0 = 0;
v0 = 10.0;
e_y0 = 0.0;
e_p0 = 0.0;

lw = 2.0; % Wheelbase
steering_lim = 0.5;
c = 0.1; % Viscous friction of the longitudinal dynamics
throttle_lim = 1.0; % Longitudinal force limit
max_speed = 50.0;
min_speed = 0.1;

%% Load
load yas_marina_full.mat

%% Setting Optimization Problem
size_vec = size(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr);
N = size_vec(1) - 1;
opti = casadi.Opti(); % Optimization problem

% ---- decision variables ---------
X = opti.variable(4, N+1); % state trajectory in path frame
t = X(1,:);
v = X(2,:);
ey = X(3,:);
ep = X(4,:);

U = opti.variable(2,N);   % throttle and steering

% ---- objective          ---------
% opti.minimize(t(end)); % minimize time
opti.minimize(0.0 * dot(U(2,:), U(2,:)) + t(end)); % minimize steering

% ---- dynamic constraints --------
% x' = [t, v, ey, ep]
f = @(x,u,kappa) [
    (1 - kappa * x(3)) / (x(2) * cos(x(4)));
    (u(1) - c * x(2)) * (1 - kappa * x(3)) / (x(2) * cos(x(4)));
    (x(2) * sin(x(4))) * ((1 - kappa * x(3)) / (x(2) * cos(x(4))));
    (x(2) * tan(u(2)) / lw) * ((1 - kappa * x(3)) / (x(2) * cos(x(4)))) - kappa; 
   ];

guess_steering = [];

for k=1:N % loop over control intervals
   kappa = x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k,4);
   ds = x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k + 1, 5) ...
       - x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k, 5);

   guess_steering = [guess_steering, tan(kappa * lw)];
   
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

   opti.subject_to(X(3,k+1) <= x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k + 1, 10));
   opti.subject_to(X(3,k+1) >= x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(k + 1, 11));

   opti.subject_to(X(1,k+1) > X(1,k));
end

% Steering is limited
opti.subject_to(-steering_lim <= U(2,:))
opti.subject_to(U(2,:) <= steering_lim);

% Throttle limit
opti.subject_to(-throttle_lim <= U(1,:))
opti.subject_to(U(1,:) <= throttle_lim);

% Speed Limits
opti.subject_to(X(2,:) <= max_speed);
opti.subject_to(X(2,:) >= min_speed);

% Angle Limits
opti.subject_to(X(4,:) <= pi/2.2);
opti.subject_to(X(4,:) >= -pi/2.2);

% ---- Initial guess
% Problem is symmetric when the pendulum is pointing down or up
%   adding a random initial guess helps it to solve
opti.set_initial(X(2,:), v0);
opti.set_initial(U(2,:), guess_steering);

% ---- initial conditions --------
opti.subject_to(t(1) == t0);
opti.subject_to(v(1) == v0);
opti.subject_to(ey(1) == e_y0);
opti.subject_to(ep(1) == e_p0);

% ---- final conditions -------
opti.subject_to(v(end) == v0);
opti.subject_to(ey(end) == e_y0);
opti.subject_to(ep(end) == e_p0);

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
opts.ipopt.max_iter = 10000;
opti.solver('ipopt', opts); % set numerical backend
sol = opti.solve();   % actual solve

%% Extract Input Sequence and Save
distance_arr = (0:N) * ds;
time_arr = sol.value(t);
vel_arr = sol.value(v);
lat_err_arr = sol.value(ey);
ang_err_arr = sol.value(ep);

distance_input_arr = (0:N-1) * ds;
steering_arr = sol.value(U(2,:));
throttle_arr = sol.value(U(1,:));

steering_ts = timeseries(steering_arr, time_arr(1:end-1));
throttle_ts = timeseries(throttle_arr, time_arr(1:end-1));

%% Transform Predictions to Map
traj_x = x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,1) + lat_err_arr' ...
    .* -sin(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,3));
traj_y = x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,2) + lat_err_arr' ...
    .* cos(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,3));

%% Post Process
set(0,'DefaultFigureWindowStyle','docked')

figure(1)
subplot(2,2,1)
hold on
plot(distance_arr, lat_err_arr, 'LineWidth', 2)
plot(distance_arr, x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:, 10), '-.', 'LineWidth', 2)
plot(distance_arr, x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:, 11), '-.', 'LineWidth', 2)
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
plot(time_arr(1:end-1), sol.value(throttle_arr), 'LineWidth', 2);
legend('$\delta$', 'Throttle', 'Location', 'northeast', 'Interpreter', 'latex')
title('Time ve Inputs')
grid minor

ts = 0.01;
out = sim("lateral_pendulum_planar_moving_vehicle_sim");
load yas_marina_min_curve
figure(3)
hold on
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,1), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,2), 'LineWidth', 2);
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,6), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,7), 'LineWidth', 2);
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,8), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,9), 'LineWidth', 2);
plot(traj_x, traj_y, 'LineWidth', 2);
plot(out.simout(:,1), out.simout(:,2), 'LineWidth', 2);
plot(yas_marina_min_curve(:,1),yas_marina_min_curve(:,2), 'LineWidth', 2);
title('Map')
legend('CenterLine','LeftMargin','RightMargin','OptimalTrajectory',...
    'OptimalInputsOpenLoopSim', 'MinCurvature')
grid minor
title('Map')
daspect([1,1,1])