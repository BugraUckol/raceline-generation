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

X = [u,v,w,p,q,r,s,en,eb,ephi,ethe,epsi]

U = [Fz, Mx, My, Mz]
        ________
-[X]-->|        |
       | System |---[X]->
-[U]-->|________|

%}

%% Constants
Fz_lim = 10;
M_lim = 10.0;
I_xx = 1;
I_yy = 1;
I_zz = 1;
m = 1;

%% Initial conditions
u0 = 1.0;
v0 = 0.0;
w0 = 0.0;

p0 = 0.0;
q0 = 0.0;
r0 = 0.0;

t0 = 0;
e_n0 = 0.1;
e_b0 = 0.0;

e_phi0 = 0.0;
e_psi0 = 0.0;
e_the0 = 0.0;

%% Setting Optimization Problem
size_vec = floor(size(path.s_arr));
N = 150; % size_vec(2) - 1;

opti = casadi.Opti();

X = opti.variable(12, N+1); % state trajectory in path frame
t = X(1,:);
u = X(2,:);
v = X(3,:);
w = X(4,:);
p = X(5,:);
q = X(6,:);
r = X(7,:);
en = X(8,:);
eb = X(9,:);
ephi = X(10,:);
ethe = X(11,:);
epsi = X(12,:);

U = opti.variable(4,N);   % Angular rates
Fz = U(1,:);
Mx = U(2,:);
My = U(3,:);
Mz = U(4,:);

% Cost Function
% Minimize angular velocity inputs p and q
% opti.minimize(1.0 * U(1,:) * U(1,:)' + 1.0 * U(2,:) * U(2,:)');
opti.minimize(10.0 * (en * en') + 1.0 * (eb * eb') + ...
    0.1 * (p * p') + 0.1 * (q * q') + 0.1 * (r * r'));
% opti.minimize(1.0 * (t * t'));

% System dynamics
% x* = [t,u,v,w,p,q,r,en,eb,ephi,ethe,epsi] states in spatial formulation
f = @(tt,u,v,w,p,q,r,en,eb,ephi,ethe,epsi,Fz,Mx,My,Mz,kappa,tau) [
    -(en*kappa - 1)/(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe));
    ((en*kappa - 1)*(q*w - r*v))/(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe));
    -((en*kappa - 1)*(p*w - r*u))/(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe));
    -((en*kappa - 1)*(q*u - p*v + Fz/m))/(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe));
    -((en*kappa - 1)*(Mx + I_yy*q*r - I_zz*q*r))/(I_xx*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)));
    -((en*kappa - 1)*(My - I_xx*p*r + I_zz*p*r))/(I_yy*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)));
    -((en*kappa - 1)*(Mz + I_xx*p*q - I_yy*p*q))/(I_zz*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)));
    -((en*kappa - 1)*(v*(cos(ephi)*cos(epsi) + sin(ephi)*sin(epsi)*sin(ethe)) - w*(cos(epsi)*sin(ephi) - cos(ephi)*sin(epsi)*sin(ethe)) + u*cos(ethe)*sin(epsi) - (eb*tau*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)))/(en*kappa - 1)))/(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe));
    -((en*kappa - 1)*(w*cos(ephi)*cos(ethe) - u*sin(ethe) + v*cos(ethe)*sin(ephi) + (en*tau*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)))/(en*kappa - 1)))/(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe));
    -((en*kappa - 1)*(p*cos(ethe) + r*cos(ephi)*sin(ethe) + q*sin(ephi)*sin(ethe) + (tau*cos(epsi)*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)))/(en*kappa - 1)))/(cos(ethe)*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)));
    ((en*kappa - 1)*(r*sin(ephi) - q*cos(ephi) + (tau*sin(epsi)*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)))/(en*kappa - 1)))/(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe));
    -((en*kappa - 1)*(r*cos(ephi) + q*sin(ephi) + (kappa*cos(ethe)*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)))/(en*kappa - 1) + (tau*cos(epsi)*sin(ethe)*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)))/(en*kappa - 1)))/(cos(ethe)*(w*(sin(ephi)*sin(epsi) + cos(ephi)*cos(epsi)*sin(ethe)) - v*(cos(ephi)*sin(epsi) - cos(epsi)*sin(ephi)*sin(ethe)) + u*cos(epsi)*cos(ethe)));
    ];

for k=1:N % loop over control intervals
   kappa = path.kappa_arr(k);
   tau = path.tau_arr(k);
   ds = path.s_arr(k + 1) ...
       - path.s_arr(k);
   
   % 1st Order Explicit Euler's Integration
   k1 = f(X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), X(6,k),...
       X(7,k), X(8,k), X(9,k), X(10,k), X(11,k), X(12,k),...
       U(1,k), U(2,k), U(3,k), U(4,k), kappa, tau);
   x_next = X(:,k) + ds * k1;
   
   % Multiple shooting
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end
opti.subject_to(t(2:N+1) > t(1:N)) % Time must increase!

% Constraints
opti.subject_to(U(1,:) <= Fz_lim);
opti.subject_to(U(2,:) <= M_lim);
opti.subject_to(U(3,:) <= M_lim);
opti.subject_to(U(4,:) <= M_lim);
opti.subject_to(U(1,:) >= -Fz_lim);
opti.subject_to(U(2,:) >= -M_lim);
opti.subject_to(U(3,:) >= -M_lim);
opti.subject_to(U(4,:) >= -M_lim);

opti.subject_to(X(8,:) <= 1.2);
opti.subject_to(X(8,:) >= -1.2);
opti.subject_to(X(9,:) <= 1.2);
opti.subject_to(X(9,:) >= -1.2);

opti.subject_to(t(1) == 0.0);
opti.subject_to(u(1) == u0);
opti.subject_to(v(1) == v0);
opti.subject_to(w(1) == w0);
opti.subject_to(p(1) == p0);
opti.subject_to(q(1) == q0);
opti.subject_to(r(1) == r0);
opti.subject_to(en(1) == e_n0);
opti.subject_to(eb(1) == e_b0);
opti.subject_to(ephi(1) == e_phi0);
opti.subject_to(ethe(1) == e_the0);
opti.subject_to(epsi(1) == e_psi0);

opti.subject_to(en(end) == 0.0);
opti.subject_to(eb(end) == 0.0);

% Initial Guess
opti.set_initial(u, 1);

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
u_arr = sol.value(u);
v_arr = sol.value(v);
w_arr = sol.value(w);
p_arr = sol.value(p);
q_arr = sol.value(q);
r_arr = sol.value(r);
en_arr = sol.value(en);
eb_arr = sol.value(eb);
ephi_arr = sol.value(ephi);
ethe_arr = sol.value(ethe);
epsi_arr = sol.value(epsi);

Fz_arr = sol.value(Fz);
Mx_arr = sol.value(Mx);
My_arr = sol.value(My);
Mz_arr = sol.value(Mz);

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
subplot(2,3,[4,5])
plot(distance_arr, ephi_arr), hold on
plot(distance_arr, ethe_arr)
plot(distance_arr, epsi_arr)
title('dist vs angles'), legend('phi', 'the', 'psi')
subplot(2,3,6)
plot(distance_arr(1:end-1), Fz_arr); hold on;
plot(distance_arr(1:end-1), Mx_arr);
plot(distance_arr(1:end-1), My_arr);
plot(distance_arr(1:end-1), Mz_arr);
legend('Fz', 'Mx', 'My', 'Mz')
title('commands')

figure(2)
plot3(path.x_arr, path.y_arr, path.z_arr, 'LineWidth', 2); hold on
plot3(x_arr, y_arr, z_arr, 'LineWidth', 2); grid minor
daspect([1,1,1])
p0 = [path.x_arr(1); path.y_arr(1); path.z_arr(1)];
dcm_b_p = CB2E([e_phi0, e_the0, e_psi0]);
dcm_p_e = CB2E([path.roll_arr(1), path.pitch_arr(1), path.yaw_arr(1)]);
r = dcm_p_e * [0; e_n0; e_b0] + p0;
% dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
dcm_b_e  = dcm_p_e * dcm_b_p;
r_hist = [r];
edges = [0, 0, 0, 0, 0;
         1.2, -1.2, -1.2, 1.2, 1.2;
         1.2, 1.2, -1.2, -1.2, 1.2];
edges_e_list = [];
for k=1:1:N-1
    cla;
    plot3(path.x_arr, path.y_arr, path.z_arr, 'LineWidth', 2)
    plot3(x_arr, y_arr, z_arr, 'LineWidth', 2)
    % plot3(r_hist(1,:), r_hist(2,:), r_hist(3,:), 'LineWidth', 3)
    
    px = [r, r + dcm_b_e(:,1)];
    py = [r, r + dcm_b_e(:,2)];
    pz = [r, r + dcm_b_e(:,3)];
    
    % plot3(px(1,:), px(2,:), px(3,:), 'LineWidth', 2, 'Color', 'r');
    % plot3(py(1,:), py(2,:), py(3,:), 'LineWidth', 2, 'Color', 'g');
    % plot3(pz(1,:), pz(2,:), pz(3,:), 'LineWidth', 2, 'Color', 'b');

    dt = time_arr(k+1) - time_arr(k);
    ds = distance_arr(k+1) - distance_arr(k);

    r = r + dcm_b_e * [u_arr(k); v_arr(k); w_arr(k)] * dt;

    r_hist = [r_hist, r];

    w_be = [p_arr(k), q_arr(k), r_arr(k)];
    dcm_b_e = dcm_b_e + dt * dcm_b_e * skew(w_be);
    dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
    pt = [path.x_arr(k); path.y_arr(k); path.z_arr(k)];
    ptx = [pt, pt + dcm_p_e(:,1) * 1.5];
    pty = [pt, pt + dcm_p_e(:,2) * 1.5];
    ptz = [pt, pt + dcm_p_e(:,3) * 1.5];
    plot3(ptx(1,:), ptx(2,:), ptx(3,:), 'LineWidth', 3, 'Color', 'r');
    plot3(pty(1,:), pty(2,:), pty(3,:), 'LineWidth', 3, 'Color', 'g');
    plot3(ptz(1,:), ptz(2,:), ptz(3,:), 'LineWidth', 3, 'Color', 'b');

    edges_e = pt + dcm_p_e * edges;
    % load edges_e_list
    % for ie = 1:5:998
    %     edges_ei = edges_e_list(:,:,ie);
    %     plot3(edges_ei(1,:), edges_ei(2,:), edges_ei(3,:), 'k', 'LineWidth', 0.15)
    % end
    % plot3(squeeze(edges_e_list(1,1,:)), squeeze(edges_e_list(2,1,:)), squeeze(edges_e_list(3,1,:)), 'k', 'LineWidth', 1.0)
    % plot3(squeeze(edges_e_list(1,2,:)), squeeze(edges_e_list(2,2,:)), squeeze(edges_e_list(3,2,:)), 'k', 'LineWidth', 1.0)
    % plot3(squeeze(edges_e_list(1,3,:)), squeeze(edges_e_list(2,3,:)), squeeze(edges_e_list(3,3,:)), 'k', 'LineWidth', 1.0)
    % plot3(squeeze(edges_e_list(1,4,:)), squeeze(edges_e_list(2,4,:)), squeeze(edges_e_list(3,4,:)), 'k', 'LineWidth', 1.0)
    % plot3(edges_e(1,:), edges_e(2,:), edges_e(3,:), 'k', 'LineWidth', 0.5)
    dcm_gt_p = CB2E([ephi_arr(k), ethe_arr(k), epsi_arr(k)]);
    dcm_gt_e = dcm_p_e * dcm_gt_p;
    rgt = [x_arr(k); y_arr(k); z_arr(k)];
    pgtx = [rgt, rgt + dcm_gt_e(:,1) * 1.5];
    pgty = [rgt, rgt + dcm_gt_e(:,2) * 1.5];
    pgtz = [rgt, rgt + dcm_gt_e(:,3) * 1.5];
    
    plot3(pgtx(1,:), pgtx(2,:), pgtx(3,:), 'LineWidth', 3, 'Color', 'r');
    plot3(pgty(1,:), pgty(2,:), pgty(3,:), 'LineWidth', 3, 'Color', 'g');
    plot3(pgtz(1,:), pgtz(2,:), pgtz(3,:), 'LineWidth', 3, 'Color', 'b');
    drawnow;
    daspect([1,1,1])
    pause(0.05)
end

% plot3(r_hist(1,:), r_hist(2,:), r_hist(3,:), 'LineWidth', 3)