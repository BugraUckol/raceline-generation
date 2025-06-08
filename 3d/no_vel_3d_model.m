%% Prep
clc, clear, close all
% Clear is especially important since the opti object should not be 
%   given constraints from the past

%% Casadi Imports
% addpath("/Users/bugrauckol/Documents/share/casadi-3")
addpath("C:\Program Files\casadi-3.6.7-windows64-matlab2018b")
import casadi.*

%% Import path properties
path = load('three_d_infinity.mat');
% path = load('circle_2d.mat');

%% System Model
%{

X = [t, ey, ez, e_psi, e_the, e_phi]
    - Time
    - Error in y direction of the Frenet Frame
    - Error in z direction of the Frenet Frame
    - Yaw of the Body Frame wrt. Frenet Frame
    - Pitch of the Body Frame wrt. Frenet Frame
    - Roll of the Body Frame wrt. Frenet Frame

U = [p, q, r] (all wrt. Earth)
    - Roll rate
    - Pitch rate
    - Yaw rate
        ________
-[X]-->|        |
       | System |---[X]->
-[U]-->|________|

%}

%% Constants
v0 = 1.0; % Not a state for constant velocity model
pqr_lim = 0.8;

%% Initial conditions
t0 = 0;
e_n0 = 0.0;
e_b0 = 0.0;
e_phi0 = 0.0;
e_psi0 = 0.0;
e_the0 = 0.0;

%% Setting Optimization Problem
size_vec = floor(size(path.s_arr));
N = size_vec(2) - 1;

opti = casadi.Opti();

X = opti.variable(6, N+1); % state trajectory in path frame
t = X(1,:);
en = X(2,:);
eb = X(3,:);
ephi = X(4,:);
ethe = X(5,:);
epsi = X(6,:);

U = opti.variable(3,N);   % Angular rates
p_com = U(1,:);
q_com = U(2,:);
r_com = U(3,:);

% Cost Function
% Minimize angular velocity inputs p and q
% opti.minimize(1.0 * U(1,:) * U(1,:)' + 1.0 * U(2,:) * U(2,:)');
% opti.minimize(1.0 * (en * en') + 1.0 * (eb * eb'));

opti.minimize((t * t'));
opti.subject_to(X(2,:) <= 2.0);
opti.subject_to(X(2,:) >= -2.0);
opti.subject_to(X(3,:) <= 2.0);
opti.subject_to(X(3,:) >= -2.0);
opti.subject_to(X(2,end) == 0.0);
opti.subject_to(X(3,end) == 0.0);

% System dynamics
% x* = [t, ey, ep, e_phi, e_the, e_psi] states in spatial formulation
% ERROR. Negating the first term of the eeb eqution seems like the fix
f = @(tt,een,eeb,eephi,eethe,eepsi,p,q,r,kappa,tau) [
    (1 - kappa * een) / (v0 * cos(eepsi) * cos(eethe));
    (1 - kappa * een) * tan(eepsi) + tau * eeb;
    (1 - kappa * een) * tan(eethe) / cos(eepsi) - tau * een; %% NEGATE
    -(tau*v0*cos(eepsi)^2*cos(eethe) - r*cos(eephi)*sin(eethe) - q*sin(eephi)*sin(eethe) - p*cos(eethe) + een*kappa*p*cos(eethe) + een*kappa*r*cos(eephi)*sin(eethe) + een*kappa*q*sin(eephi)*sin(eethe))/(v0*cos(eepsi)*cos(eethe)^2);
    (q*cos(eephi) - r*sin(eephi) - een*kappa*q*cos(eephi) + een*kappa*r*sin(eephi) + tau*v0*cos(eepsi)*cos(eethe)*sin(eepsi))/(v0*cos(eepsi)*cos(eethe));
    -(kappa*v0*cos(eepsi)*cos(eethe)^2 - q*sin(eephi) - r*cos(eephi) + een*kappa*r*cos(eephi) + een*kappa*q*sin(eephi) + tau*v0*cos(eepsi)^2*cos(eethe)*sin(eethe))/(v0*cos(eepsi)*cos(eethe)^2)
];

for k=1:N % loop over control intervals
   kappa = path.kappa_arr(k);
   tau = path.tau_arr(k);
   ds = path.s_arr(k + 1) ...
       - path.s_arr(k);
   
   % 1st Order Explicit Euler's Integration
   k1 = f(X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), X(6,k),...
       U(1,k), U(2,k), U(3,k), kappa, tau);
   x_next = X(:,k) + ds * k1;
   
   % Multiple shooting
   opti.subject_to(X(:,k+1)==x_next); % close the gaps

   if k < 2
    opti.subject_to(U(1,k) == 0.0);
    opti.subject_to(U(2,k) == 0.0);
    opti.subject_to(U(3,k) == 0.0);
   end
end
opti.subject_to(t(2:N+1) > t(1:N)) % Time must increase!

% Constraints
opti.subject_to(U(1,:) <= pqr_lim);
opti.subject_to(U(2,:) <= pqr_lim);
opti.subject_to(U(3,:) <= pqr_lim);
opti.subject_to(U(1,:) >= -pqr_lim);
opti.subject_to(U(2,:) >= -pqr_lim);
opti.subject_to(U(3,:) >= -pqr_lim);

opti.subject_to(X(5,:) <= 1.2);
opti.subject_to(X(5,:) >= -1.2);

% opti.subject_to(U(3,:) == 0.0);

opti.subject_to(t(1) == 0.0);
opti.subject_to(en(1) == e_n0);
opti.subject_to(eb(1) == e_b0);
opti.subject_to(ephi(1) == e_phi0);
opti.subject_to(ethe(1) == e_the0);
opti.subject_to(epsi(1) == e_psi0);

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
ephi_arr = sol.value(ephi);
ethe_arr = sol.value(ethe);
epsi_arr = sol.value(epsi);

p_com_arr = sol.value(p_com);
q_com_arr = sol.value(q_com);
r_com_arr = sol.value(r_com);

%% 3D Recreation
x_arr = zeros(1,length(distance_arr));
y_arr = zeros(1,length(distance_arr));
z_arr = zeros(1,length(distance_arr));

k = 0;
for point = distance_arr
    k = k + 1;

    c_p2e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
    rpe_e = [path.x_arr(k); path.y_arr(k); path.z_arr(k)];
    rbp_e = c_p2e * [0; en_arr(k); eb_arr(k)];

    rbe_e = rpe_e + rbp_e;

    x_arr(k) = rbe_e(1,1);
    y_arr(k) = rbe_e(2,1);
    z_arr(k) = rbe_e(3,1);
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
plot(distance_arr(1:end-1), r_com_arr);
title('commands')

figure(2)
plot3(path.x_arr, path.y_arr, path.z_arr, 'LineWidth', 3); hold on
plot3(x_arr, y_arr, z_arr, 'LineWidth', 3); grid minor
daspect([1,1,1])
p0 = [path.x_arr(1); path.y_arr(1); path.z_arr(1)];
dcm_b_p = CB2E([e_phi0, e_the0, e_psi0]);
dcm_p_e = CB2E([path.roll_arr(1), path.pitch_arr(1), path.yaw_arr(1)]);
r = dcm_p_e * [0; e_n0; e_b0] + p0;
% dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
dcm_b_e  = dcm_p_e * dcm_b_p;
r_hist = [r];
for k=1:1:N-1
    cla;
    plot3(path.x_arr, path.y_arr, path.z_arr, 'LineWidth', 3)
    plot3(x_arr, y_arr, z_arr, 'LineWidth', 3)
    plot3(r_hist(1,:), r_hist(2,:), r_hist(3,:), 'LineWidth', 3)
    
    px = [r, r + dcm_b_e(:,1)];
    py = [r, r + dcm_b_e(:,2)];
    pz = [r, r + dcm_b_e(:,3)];
    
    plot3(px(1,:), px(2,:), px(3,:), 'LineWidth', 2, 'Color', 'r');
    plot3(py(1,:), py(2,:), py(3,:), 'LineWidth', 2, 'Color', 'g');
    plot3(pz(1,:), pz(2,:), pz(3,:), 'LineWidth', 2, 'Color', 'b');

    dt = time_arr(k+1) - time_arr(k);
    ds = distance_arr(k+1) - distance_arr(k);

    r = r + dcm_b_e * [v0; 0; 0] * dt;

    r_hist = [r_hist, r];

    w_be = [p_com_arr(k), q_com_arr(k), r_com_arr(k)]
    r_bp_p = [0; en_arr(k); eb_arr(k)]
    dcm_b_e = dcm_b_e + dt * dcm_b_e * skew(w_be);
    dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
    pt = [path.x_arr(k); path.y_arr(k); path.z_arr(k)];
    ptx = [pt, pt + dcm_p_e(:,1) * 0.5];
    pty = [pt, pt + dcm_p_e(:,2) * 0.5];
    ptz = [pt, pt + dcm_p_e(:,3) * 0.5];
    plot3(ptx(1,:), ptx(2,:), ptx(3,:), 'LineWidth', 2, 'Color', 'r');
    plot3(pty(1,:), pty(2,:), pty(3,:), 'LineWidth', 2, 'Color', 'g');
    plot3(ptz(1,:), ptz(2,:), ptz(3,:), 'LineWidth', 2, 'Color', 'b');
    
    dcm_gt_p = CB2E([ephi_arr(k), ethe_arr(k), epsi_arr(k)]);
    dcm_gt_e = dcm_p_e * dcm_gt_p;
    rgt = [x_arr(k); y_arr(k); z_arr(k)];
    pgtx = [rgt, rgt + dcm_gt_e(:,1) * 0.5];
    pgty = [rgt, rgt + dcm_gt_e(:,2) * 0.5];
    pgtz = [rgt, rgt + dcm_gt_e(:,3) * 0.5];
    
    plot3(pgtx(1,:), pgtx(2,:), pgtx(3,:), 'LineWidth', 2, 'Color', 'r');
    plot3(pgty(1,:), pgty(2,:), pgty(3,:), 'LineWidth', 2, 'Color', 'g');
    plot3(pgtz(1,:), pgtz(2,:), pgtz(3,:), 'LineWidth', 2, 'Color', 'b');

    drawnow;
    daspect([1,1,1])
    pause(0.05)
end

plot3(r_hist(1,:), r_hist(2,:), r_hist(3,:), 'LineWidth', 3)