%% Casadi Imports
addpath("C:\Program Files\casadi-3.6.7-windows64-matlab2018b")
import casadi.*

%% Initial conditions and constants
t0 = 0;
v0 = 1.0;
e_y0 = 0.0;
e_p0 = 0.0;

lw = 2.0; % Wheelbase
steering_lim = 0.5;
c = 0.1; % Viscous friction of the longitudinal dynamics
throttle_lim = 1.0; % Longitudinal force limit
max_speed = 10.0;
min_speed = 1.0;

%% Setting Optimization Problem
N = 5000;
opti = casadi.Opti(); % Optimization problem

% ---- decision variables ---------
X = opti.variable(4, N+1); % state trajectory in path frame
t = X(1,:);
v = X(2,:);
ey = X(3,:);
ep = X(4,:);

U = opti.variable(2,N);   % throttle and steering

% ---- objective          ---------
opti.minimize(t(end)); % hold the pendulum straight

% ---- dynamic constraints --------
% x' = [t, v, ey, ep]
f = @(x,u,kappa) [
    (1 - kappa * x(3)) / (x(2) * cos(x(4)));
    (u(1) - c * x(2)) * (1 - kappa * x(3)) / (x(2) * cos(x(4)));
    (x(2) * sin(x(4))) * ((1 - kappa * x(3)) / (x(2) * cos(x(4))));
    (x(2) * tan(u(2)) / lw) * ((1 - kappa * x(3)) / (x(2) * cos(x(4)))) - kappa; 
   ];

for k=1:N % loop over control intervals
   kappa = 0.1;
   ds = 0.01;
   % Runge-Kutta 4 integration
   k1 = f(X(:,k),         U(:,k), kappa);
   k2 = f(X(:,k)+ds/2*k1, U(:,k), kappa);
   k3 = f(X(:,k)+ds/2*k2, U(:,k), kappa);
   k4 = f(X(:,k)+ds*k3,   U(:,k), kappa);
   x_next = X(:,k) + ds/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
   % if k < N
   %  opti.subject_to((U(:,k+1) - U(:,k))/dt < 10.0);
   %  opti.subject_to((U(:,k+1) - U(:,k))/dt > -10.0);
   % end
   % k1 = f(X(:,k), U(:,k));
   % x_next = X(:,k) + dt*k1;
   % opti.subject_to(X(:,k+1)==x_next);
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
opti.subject_to(X(3,:) <= 1.0);
opti.subject_to(X(3,:) >= -1.0);

% ---- Initial guess
% Problem is symmetric when the pendulum is pointing down or up
%   adding a random initial guess helps it to solve
opti.set_initial(X(2,:), v0);

% ---- initial conditions --------
opti.subject_to(t(1) == t0);
opti.subject_to(v(1) == v0);
opti.subject_to(ey(1) == e_y0);
opti.subject_to(ep(1) == e_p0);

% ---- final conditions -------
% opti.subject_to(psi(end) == 0.0);
% opti.subject_to(pos_y(end) == 0.0);
% opti.subject_to(theta(end) == 0.0);
% opti.subject_to(theta_d(end) == 0.0);

% ---- set initial guess
% opti.set_initial(U, 0.5);

%% Solve NLP              ------
opti.solver('ipopt'); % set numerical backend
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

%% Post Process
figure(1)
subplot(2,2,1)
hold on
plot(distance_arr, lat_err_arr, 'LineWidth', 2)
title('Lat Error')
grid minor
subplot(2,2,2)
hold on
plot(distance_arr, ang_err_arr, 'LineWidth', 2)
title('Ang Error')
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
title('Time')
grid minor

figure(2)
hold on
plot(time_arr(1:end-1), sol.value(steering_arr), 'LineWidth', 2);
plot(time_arr(1:end-1), sol.value(throttle_arr), 'LineWidth', 2);
legend('$\delta$', 'Throttle', 'Location', 'northeast', 'Interpreter', 'latex')
grid minor
