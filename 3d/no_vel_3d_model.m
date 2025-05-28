%% Prep
clc, clear, close all
% Clear is especially important since the opti object should not be 
%   given constraints from the past

%% Casadi Imports
addpath("/Users/bugrauckol/Documents/share/casadi-3")
% addpath("C:\Program Files\casadi-3.6.7-windows64-matlab2018b")
import casadi.*

%% Import path properties
path = load('three_d_infinity.mat');

%% System Model
%{

X = [t, ey, ez, e_psi, e_the, e_phi]
    - Time
    - Error in y direction of the Frenet Frame
    - Error in z direction of the Frenet Frame
    - Yaw of the Body Frame wrt. Frenet Frame
    - Pitch of the Body Frame wrt. Frenet Frame
    - Roll of the Body Frame wrt. Frenet Frame
        ________
-[X]-->|        |
-[p]-->| System |---[X]->
-[q]-->|________|

%}

%% Constants
v0 = 1.0; % Not a state for constant velocity model
pq_lim = 0.8;

%% Initial conditions
t0 = 0;
e_y0 = 0.0;
e_z0 = 0.0;
e_psi = 0.0;
e_the = 0.0;
e_phi = 0.0;

%% Setting Optimization Problem
size_vec = floor(size(path.s_arr));
N = size_vec(2) - 1;

opti = casadi.Opti();

X = opti.variable(5, N+1); % state trajectory in path frame
t = X(1,:);
en = X(2,:);
eb = X(3,:);
epsi = X(4,:);
ethe = X(5,:);

U = opti.variable(2,N);   %steering
p_com = U(1,:);
q_com = U(2,:);

% Cost Function
% opti.minimize(1.0 * U(1,:) * U(1,:)' + 1.0 * U(2,:) * U(2,:)');
opti.minimize(1.0 * X(2,:) * X(2,:)' + 1.0 * X(3,:) * X(3,:)' ...
    + 0.0 * X(4,:) * X(4,:)' + + 0.0 * X(5,:) * X(5,:)');
% Minimize angular velocity inputs p and q

% ---- dynamic constraints --------
% x' = [t, ey, ep, e_the, e_psi]
f = @(tt,een,eeb,eepsi,eethe,p,q,kappa,tau) [
    (1 - kappa * een) / (v0 * cos(eepsi) * cos(eethe));
    (1 - kappa * een) * tan(eepsi) + tau * eeb;
    (1 - kappa * een) * tan(eethe)/cos(eepsi) - tau * een;
    p - tau * sin(eepsi);
    q - kappa + tau*tan(eethe) - tau * sin(eepsi) * tan(eepsi) * tan(eethe) 
   ];

for k=1:N % loop over control intervals
   kappa = path.kappa_arr(k);
   tau = path.tau_arr(k);
   ds = path.s_arr(k + 1) ...
       - path.s_arr(k);
   
   % 1st Order Explicit Euler's Integration
   k1 = f(X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), U(1,k), U(2,k), kappa, tau);
   x_next = X(:,k) + ds * k1;
   
   % Multiple shooting
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

% Constraints
opti.subject_to(U(1,:) <= pq_lim);
opti.subject_to(U(2,:) <= pq_lim);
opti.subject_to(U(1,:) >= -pq_lim);
opti.subject_to(U(2,:) >= -pq_lim);

opti.subject_to(t(1) == 0.0);
opti.subject_to(en(1) == 0.3);
opti.subject_to(eb(1) == 0.0);
opti.subject_to(epsi(1) == 0.2);
opti.subject_to(ethe(1) == 0.2); %%should be nonzero!!! in this convention

%% Solve the problem
opts = struct();
opts.ipopt = struct();
opts.ipopt.tol = 1e-3;  % Increase tolerance
opts.ipopt.acceptable_tol = 1e-4;  % Reduce acceptable tolerance
opts.ipopt.max_iter = 10000;
opts.ipopt.mu_init = 1e-3;  % Initial value for the barrier parameter
opti.solver('ipopt', opts); % Set numerical backend
sol = opti.solve();   % actual solve

%% Extract Solutions
distance_arr = path.s_arr(1:N+1);
time_arr = sol.value(t);
en_arr = sol.value(en);
eb_arr = sol.value(eb);
epsi_arr = sol.value(epsi);
ethe_arr = sol.value(ethe);

p_com_arr = sol.value(p_com);
q_com_arr = sol.value(q_com);

%% 3D Recreation
x_arr = [];
y_arr = [];
z_arr = [];

k = 0;
for point = distance_arr
    k = k + 1;

    dcm = angle2dcm(path.yaw_arr(k), path.pitch_arr(k), path.roll_arr(k));
    rpe_vec = [path.x_arr(k), path.y_arr(k), path.z_arr(k)]';
    rbp_vec = dcm * [0, en_arr(k), eb_arr(k)]';

    rbe_vec = rpe_vec + rbp_vec;

    x_arr = [x_arr, rbe_vec(1,1)];
    y_arr = [y_arr, rbe_vec(2,1)];
    z_arr = [z_arr, rbe_vec(3,1)];
end

%% Plots
set(0,'DefaultFigureWindowStyle','docked')
figure(1)
subplot(2,3,1)
plot(distance_arr, time_arr);
title('dist vs time')
subplot(2,3,2)
plot(distance_arr, en_arr)
title('dist vs en')
subplot(2,3,3)
plot(distance_arr, eb_arr)
title('dist vs eb')
subplot(2,3,4)
plot(distance_arr, epsi_arr)
title('dist vs epsi')
subplot(2,3,5)
plot(distance_arr, ethe_arr)
title('dist vs ethe')
subplot(2,3,6)
plot(distance_arr(1:end-1), p_com_arr); hold on;
plot(distance_arr(1:end-1), q_com_arr);
title('commands')

figure(2)
plot3(path.x_arr, path.y_arr, path.z_arr); hold on
plot3(x_arr, y_arr, z_arr); grid minor