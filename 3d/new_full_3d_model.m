%% Prep
clc, clear
% Clear is especially important since the opti object should not be 
%   given constraints from the past

%% Casadi Imports
% addpath("C:\Users\bugrauckol\Desktop\bugra\casadi")
% addpath("C:\Program Files\casadi-3.6.7-windows64-matlab2018b")
addpath("/Users/bugrauckol/Documents/share/casadi-3")
import casadi.*

%% Import path properties
% path = load('three_d_infinity.mat');
path = load('circle_2d.mat');
% path = load('trefoil.mat');
% path = load('half_circle_2d.mat');

%% System Model
%{

X = [t,en,eb,ephi,ethe,epsi,u,v,w,p,q,r]

U = [Fz, Mx, My, Mz]
        ________
-[X]-->|        |
       | System |---[X]->
-[U]-->|________|

%}

%% Constants
T_max = 5;
T_min = 0;

m = 0.250;
I_xx = 1e-2;
I_yy = I_xx;
I_zz = I_xx;
I_b = [I_xx, 0, 0; 0, I_yy, 0; 0, 0, I_zz];

Mxy_max = (T_max/4) * (0.15 * sqrt(2) / 2) * 4;
Mz_max = Mxy_max/ 5;

V_max = 10;

%% Setting Optimization Problem
size_vec = floor(size(path.s_arr));
N = size_vec(2) - 1;

opti = casadi.Opti();

X = opti.variable(12, N+1); % state trajectory in path frame
t = X(1,:);
en = X(2,:);
eb = X(3,:);
ephi = X(4,:);
ethe = X(5,:);
epsi = X(6,:);
u = X(7,:);
v = X(8,:);
w = X(9,:);
p = X(10,:);
q = X(11,:);
r = X(12,:);

U = opti.variable(4,N);   % Angular rates
T_com = U(1,:);
Mx_com = U(2,:);
My_com = U(3,:);
Mz_com = U(4,:);

% Cost Function
opti.minimize(t(end));

% System dynamics
f = @(tt,en,eb,e_phi,e_the,e_psi,u,v,w,p,q,r,Fz,Mx,My,Mz,kappa,tau) [
    % dt / ds
    -(en*kappa - 1)/(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the))
    % den / ds
    -((en*kappa - 1)*(v*(cos(e_phi)*cos(e_psi) + sin(e_phi)*sin(e_psi)*sin(e_the)) - w*(cos(e_psi)*sin(e_phi) - cos(e_phi)*sin(e_psi)*sin(e_the)) + u*cos(e_the)*sin(e_psi) - (eb*tau*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))/(en*kappa - 1)))/(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the))
    % deb / ds
    -((en*kappa - 1)*(w*cos(e_phi)*cos(e_the) - u*sin(e_the) + v*cos(e_the)*sin(e_phi) + (en*tau*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))/(en*kappa - 1)))/(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the))
    % dephi / ds (phi of the of the body frame wrt. the path frame)
    -(tau*u*cos(e_psi)^2*cos(e_the) - r*cos(e_phi)*sin(e_the) - q*sin(e_phi)*sin(e_the) - p*cos(e_the) + en*kappa*p*cos(e_the) + tau*w*cos(e_phi)*cos(e_psi)^2*sin(e_the) + tau*v*cos(e_psi)^2*sin(e_phi)*sin(e_the) + en*kappa*r*cos(e_phi)*sin(e_the) + en*kappa*q*sin(e_phi)*sin(e_the) - tau*v*cos(e_phi)*cos(e_psi)*sin(e_psi) + tau*w*cos(e_psi)*sin(e_phi)*sin(e_psi))/(cos(e_the)*(u*cos(e_psi)*cos(e_the) - v*cos(e_phi)*sin(e_psi) + w*sin(e_phi)*sin(e_psi) + w*cos(e_phi)*cos(e_psi)*sin(e_the) + v*cos(e_psi)*sin(e_phi)*sin(e_the)))
    % dethe / ds (theta of the of the body frame wrt. the path frame)
    (q*cos(e_phi) - r*sin(e_phi) - tau*v*cos(e_phi) + tau*w*sin(e_phi) + tau*v*cos(e_phi)*cos(e_psi)^2 - tau*w*cos(e_psi)^2*sin(e_phi) - en*kappa*q*cos(e_phi) + en*kappa*r*sin(e_phi) + tau*u*cos(e_psi)*cos(e_the)*sin(e_psi) + tau*w*cos(e_phi)*cos(e_psi)*sin(e_psi)*sin(e_the) + tau*v*cos(e_psi)*sin(e_phi)*sin(e_psi)*sin(e_the))/(u*cos(e_psi)*cos(e_the) - v*cos(e_phi)*sin(e_psi) + w*sin(e_phi)*sin(e_psi) + w*cos(e_phi)*cos(e_psi)*sin(e_the) + v*cos(e_psi)*sin(e_phi)*sin(e_the))
    % depsi / ds (psi of the of the body frame wrt. the path frame)
    -((en*kappa - 1)*((cos(e_phi)*(r + (tau*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the))*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))/(en*kappa - 1) + (kappa*cos(e_phi)*cos(e_the)*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))/(en*kappa - 1)))/cos(e_the) + (sin(e_phi)*(q - (tau*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the))*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))/(en*kappa - 1) + (kappa*cos(e_the)*sin(e_phi)*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))/(en*kappa - 1)))/cos(e_the)))/(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the))
    % du / ds
    ((en*kappa - 1)*(q*w - r*v))/(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the));
    % dv / ds
    -((en*kappa - 1)*(p*w - r*u))/(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the));
    % dw / ds
    -((en*kappa - 1)*(q*u - p*v + Fz/m))/(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the));
    % dp / ds
    -((en*kappa - 1)*(Mx + I_yy*q*r - I_zz*q*r))/(I_xx*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))
    % dq / ds
    -((en*kappa - 1)*(My - I_xx*p*r + I_zz*p*r))/(I_yy*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))
    % dr / ds
    -((en*kappa - 1)*(Mz + I_xx*p*q - I_yy*p*q))/(I_zz*(w*(sin(e_phi)*sin(e_psi) + cos(e_phi)*cos(e_psi)*sin(e_the)) - v*(cos(e_phi)*sin(e_psi) - cos(e_psi)*sin(e_phi)*sin(e_the)) + u*cos(e_psi)*cos(e_the)))
    ];

for k=1:N % loop over control intervals
   kappa = path.kappa_arr(k);
   tau = path.tau_arr(k);
   ds = path.s_arr(k + 1) ...
       - path.s_arr(k);

   % 1st Order Explicit Euler's Integration
   k1 = f(X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), X(6,k),...
       X(7,k), X(8,k), X(9,k), X(10,k), X(11,k), X(12,k), U(1,k), 0, 0, 0, ...
       kappa, tau);
   x_next = X(:,k) + ds * k1;
   
   % Algebraic relation between optimization parameters and propagation
   opti.subject_to(X(:,k+1)==x_next);

end

% Fundamental Constraint
opti.subject_to(t(2:N+1) > t(1:N))

% Input Constraints
opti.subject_to(T_com <= T_max);
opti.subject_to(Mx_com <= Mxy_max);
opti.subject_to(My_com <= Mxy_max);
opti.subject_to(Mz_com <= Mz_max);
opti.subject_to(T_com >= T_min);
opti.subject_to(Mx_com >= -Mxy_max);
opti.subject_to(My_com >= -Mxy_max);
opti.subject_to(Mz_com >= -Mz_max);

% Velocity Constraints
opti.subject_to(sqrt(u.^2 + v.^2 + w.^2) <= V_max);
opti.subject_to(sqrt(u.^2 + v.^2 + w.^2) >= -V_max);

% Corridor Constraints
opti.subject_to(en <= 0.5);
opti.subject_to(eb <= 0.5);
opti.subject_to(en >= -0.5);
opti.subject_to(eb >= -0.5);

% Angle Constraints
opti.subject_to(ephi <= pi);
opti.subject_to(ethe <= pi);
opti.subject_to(epsi <= pi);
opti.subject_to(ephi >= -pi);
opti.subject_to(ethe >= -pi);
opti.subject_to(epsi >= -pi);

% Initial Conditions
opti.subject_to(t(1) == 0);
% opti.subject_to(en(1) == e_n0);
% opti.subject_to(eb(1) == e_b0);
% opti.subject_to(ephi(1) == e_phi0);
% opti.subject_to(ethe(1) == e_the0);
% opti.subject_to(epsi(1) == e_psi0);
% opti.subject_to(u(1) == 0);
% opti.subject_to(v(1) == 0);
% opti.subject_to(w(1) == 0);
% opti.subject_to(p(1) == 0);
% opti.subject_to(q(1) == 0);
% opti.subject_to(r(1) == 0);

prevsol = load('prevsol_circle_uvw');
opti.set_initial(t, prevsol.time_arr);
opti.set_initial(en, prevsol.en_arr);
opti.set_initial(eb, prevsol.eb_arr);
opti.set_initial(ephi, prevsol.ephi_arr);
opti.set_initial(ethe, prevsol.ethe_arr);
opti.set_initial(epsi, prevsol.epsi_arr);
opti.set_initial(u, prevsol.u_arr);
opti.set_initial(v, prevsol.v_arr);
opti.set_initial(w, prevsol.w_arr);
opti.set_initial(p, prevsol.p_arr);
opti.set_initial(q, prevsol.q_arr);
opti.set_initial(r, prevsol.r_arr);

opti.set_initial(T_com, prevsol.T_com_arr);
opti.set_initial(Mx_com, prevsol.Mx_com_arr);
opti.set_initial(My_com, prevsol.My_com_arr);
opti.set_initial(Mz_com, prevsol.Mz_com_arr);
% opti.set_initial(u, 1);
% opti.set_initial(v, 0);
% opti.set_initial(w, 0);

% opti.subject_to(en(end) == 0.0);
% opti.subject_to(eb(end) == 0.0);

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
u_arr = sol.value(u);
v_arr = sol.value(v);
w_arr = sol.value(w);
p_arr = sol.value(p);
q_arr = sol.value(q);
r_arr = sol.value(r);

T_com_arr = sol.value(T_com);
Mx_com_arr = sol.value(Mx_com);
My_com_arr = sol.value(My_com);
Mz_com_arr = sol.value(Mz_com);

% save prevsol time_arr en_arr eb_arr ephi_arr ethe_arr epsi_arr v_arr alpha_arr beta_arr T_com_arr p_com_arr q_com_arr r_com_arr

%% 3D Recreation
x_arr = zeros(1,length(distance_arr));
y_arr = zeros(1,length(distance_arr));
z_arr = zeros(1,length(distance_arr));

k = 0;

edges = [0, 0, 0, 0, 0;
         0.5, -0.5, -0.5, 0.5, 0.5;
         0.5, 0.5, -0.5, -0.5, 0.5];
edges_e_list = [];
pp0 = [path.x_arr(1); path.y_arr(1); path.z_arr(1)];
dcm_p_e = CB2E([path.roll_arr(1), path.pitch_arr(1), path.yaw_arr(1)]);

p = dcm_p_e * [0; en_arr(1); en_arr(1)] + pp0;
% dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
for point = distance_arr
    k = k + 1;

    c_p2e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
    rpe_e = [path.x_arr(k); path.y_arr(k); path.z_arr(k)];
    rbp_e = c_p2e * [0; en_arr(k); eb_arr(k)];

    rbe_e = rpe_e + rbp_e;

    x_arr(k) = rbe_e(1,1);
    y_arr(k) = rbe_e(2,1);
    z_arr(k) = rbe_e(3,1);

    pt = [path.x_arr(k); path.y_arr(k); path.z_arr(k)];
    ptx = [pt, pt + c_p2e(:,1) * 1.5];
    pty = [pt, pt + c_p2e(:,2) * 1.5];
    ptz = [pt, pt + c_p2e(:,3) * 1.5];
    % plot3(ptx(1,:), ptx(2,:), ptx(3,:), 'LineWidth', 3, 'Color', 'r');
    % plot3(pty(1,:), pty(2,:), pty(3,:), 'LineWidth', 3, 'Color', 'g');
    % plot3(ptz(1,:), ptz(2,:), ptz(3,:), 'LineWidth', 3, 'Color', 'b');

    edges_e = pt + c_p2e * edges;
    edges_e_list = cat(3, edges_e_list, edges_e);
end

%% Plots
setDefaultFigureProperties()
set(0,'DefaultFigureWindowStyle','docked')
figure(1)
subplot(3,2,1)
plot(distance_arr, time_arr);
title('dist vs time')
hold on
subplot(3,2,2)
plot(distance_arr, en_arr), hold on
plot(distance_arr, eb_arr)
title('dist vs deviations')
legend('en', 'eb')
subplot(3,2,3)
plot(distance_arr, ephi_arr), hold on
plot(distance_arr, ethe_arr)
plot(distance_arr, epsi_arr)
title('dist vs angles'), legend('phi', 'the', 'psi')
subplot(3,2,4)
plot(distance_arr, u_arr), hold on
plot(distance_arr, v_arr) 
plot(distance_arr, w_arr) 
title('dist vs uvw')
subplot(3,2,5)
plot(distance_arr, p_arr), hold on
plot(distance_arr, q_arr)
plot(distance_arr, r_arr)
title('dist vs pqr')

figure(2)
subplot(2,1,1)
plot(distance_arr(1:end-1), T_com_arr); hold on;
title('distance vs thrust')
subplot(2,1,2)
plot(distance_arr(1:end-1), Mx_com_arr); hold on;
plot(distance_arr(1:end-1), My_com_arr);
plot(distance_arr(1:end-1), Mz_com_arr);
title('distance vs commands')
legend('p','q','r')
%%
figure(3)
subplot(2,2,1)
plot3(path.x_arr, path.y_arr, path.z_arr, 'LineWidth', 2); hold on
plot3(x_arr, y_arr, z_arr, 'LineWidth', 2); grid minor
daspect([1,1,1])
% dcm_b_p = CB2E([ephi_arr(1), ephi_arr(2), ephi_arr(3)]);
% dcm_p_e = CB2E([path.roll_arr(1), path.pitch_arr(1), path.yaw_arr(1)]);
% p = dcm_p_e * [0; en_arr(1); eb_arr(1)] + pp0;
% dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
% dcm_b_e  = dcm_p_e * dcm_b_p;
p_hist = [p];
for k=1:1:N-1
    subplot(2,2,1)
    cla;
    dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
    dcm_b_p = CB2E([ephi_arr(k), ethe_arr(k), epsi_arr(k)]);
    dcm_b_e = dcm_p_e * dcm_b_p;

    plot3(path.x_arr, path.y_arr, path.z_arr, 'LineWidth', 2)
    plot3(x_arr, y_arr, z_arr, 'LineWidth', 2)
    pp = [path.x_arr(k); path.y_arr(k); path.z_arr(k)];
    p = dcm_p_e * [0; en_arr(k); eb_arr(k)] + pp;
    
    bx = [p, p + dcm_b_e(:,1)];
    by = [p, p + dcm_b_e(:,2)];
    bz = [p, p + dcm_b_e(:,3)];
    
    plot3(bx(1,:), bx(2,:), bx(3,:), 'LineWidth', 2, 'Color', 'r');
    plot3(by(1,:), by(2,:), by(3,:), 'LineWidth', 2, 'Color', 'g');
    plot3(bz(1,:), bz(2,:), bz(3,:), 'LineWidth', 2, 'Color', 'b');

    % dt = time_arr(k+1) - time_arr(k);
    % ds = distance_arr(k+1) - distance_arr(k);

    % p = p + dcm_v_e * [v0; 0; 0] * dt;

    % p_hist = [p_hist, p];

    % w_be = [p_com_arr(k), q_com_arr(k), r_com_arr(k)];

    edges_ei = edges_e_list(:,:,k);
    plot3(edges_ei(1,:), edges_ei(2,:), edges_ei(3,:), 'k', 'LineWidth', 0.15), hold on
    plot3(squeeze(edges_e_list(1,1,:)), squeeze(edges_e_list(2,1,:)), squeeze(edges_e_list(3,1,:)), 'k', 'LineWidth', 1.0)
    plot3(squeeze(edges_e_list(1,2,:)), squeeze(edges_e_list(2,2,:)), squeeze(edges_e_list(3,2,:)), 'k', 'LineWidth', 1.0)
    plot3(squeeze(edges_e_list(1,3,:)), squeeze(edges_e_list(2,3,:)), squeeze(edges_e_list(3,3,:)), 'k', 'LineWidth', 1.0)
    plot3(squeeze(edges_e_list(1,4,:)), squeeze(edges_e_list(2,4,:)), squeeze(edges_e_list(3,4,:)), 'k', 'LineWidth', 1.0)
    plot3(edges_e(1,:), edges_e(2,:), edges_e(3,:), 'k', 'LineWidth', 0.5)
    
    % Ground truth of velocity frame
    v_vec = [u_arr(k); v_arr(k); w_arr(k)];
    vx = [p, p + dcm_b_e * v_vec / norm(v_vec)];
    % Plot Velocity Frame
    plot3(vx(1,:), vx(2,:), vx(3,:), 'LineWidth', 3, 'Color', 'm');

    % Plot Thrust Vector
    % t_vec = [p, p + dcm_b_e(:,3) * T_com_arr(k)];
    % plot3(t_vec(1,:), t_vec(2,:), t_vec(3,:), 'LineWidth', 3, 'Color', 'cyan');

    daspect([1,1,1])
    subplot(2,2,2)
    cla;
    vb = [zeros(3,1), v_vec];
    plot3(vb(1,:), vb(2,:), vb(3,:), 'LineWidth', 3, 'Color', 'm'); hold on

    bxb = [zeros(3,1), [1;0;0]];
    byb = [zeros(3,1), [0;1;0]];
    bzb = [zeros(3,1), [0;0;1]];

    plot3(bxb(1,:), bxb(2,:), bxb(3,:), 'LineWidth', 2, 'Color', 'r');
    plot3(byb(1,:), byb(2,:), byb(3,:), 'LineWidth', 2, 'Color', 'g');
    plot3(bzb(1,:), bzb(2,:), bzb(3,:), 'LineWidth', 2, 'Color', 'b');

    xlim([-10,10]);
    ylim([-10,10]);
    zlim([-10,10]);
    daspect([1,1,1])
    grid on

    subplot(2,2,3)
    cla;
    plot(distance_arr(1:end-1), T_com_arr); hold on;
    xline(distance_arr(k))
    title('distance vs thrust')
    subplot(2,2,4)
    cla;
    plot(distance_arr(1:end-1), Mx_com_arr); hold on;
    plot(distance_arr(1:end-1), My_com_arr);
    plot(distance_arr(1:end-1), Mz_com_arr);
    xline(distance_arr(k))
    title('distance vs commands')
    legend('p','q','r')

    drawnow;
    pause(0.05)

    if mod(k,10) == 0
        k;
    end
end

%%
dcm_p_e_res = CB2E([path.roll_arr(1), path.pitch_arr(1), path.yaw_arr(1)]);
dcm_b_p_res = CB2E([ephi_arr(1), ethe_arr(1), epsi_arr(1)]);
dcm_b_e_res = dcm_p_e_res * dcm_b_p_res;
pp = [path.x_arr(1); path.y_arr(1); path.z_arr(1)];
P_res = dcm_p_e_res * [0; en_arr(1); eb_arr(1)] + pp;
v_res = [u_arr(1); v_arr(1); w_arr(1)];
w_res = [p_arr(1); q_arr(1); r_arr(1)];
[a,b,c] = dcm2angle(dcm_b_e_res');
E_res = [c,b,a];
figure(3)
subplot(2,2,1)
for k=1:1:N-1
dt = time_arr(k+1) - time_arr(k);
T_res = [0;0;T_com_arr(k)];
M_res = [Mx_com_arr(k); My_com_arr(k); Mz_com_arr(k)];
C_res = CB2E(E_res);

P_res = P_res + dt * CB2E(E_res) * v_res;
v_res = v_res + dt * (T_res/m - cross(w_res, v_res));
w_res = w_res + dt * I_b^-1 * (M_res - cross(w_res, I_b * w_res));
E_res = E_res + dt * W2ED(w_res, E_res);
plot3(P_res(1), P_res(2), P_res(3), '*', 'MarkerSize', 5), hold on
end

% plot3(r_hist(1,:), r_hist(2,:), r_hist(3,:), 'LineWidth', 3)
% 
% h_ref   = plot3(NaN, NaN, NaN, 'k.-', 'LineWidth', 6);
% h_start = plot3(NaN, NaN, NaN, 'g-',  'LineWidth', 4);
% h_end   = plot3(NaN, NaN, NaN, 'r-',  'LineWidth', 4);
% h_bound = plot3(NaN, NaN, NaN, 'k-',  'LineWidth', 0.15);
% 
% legend([h_ref, h_start, h_end, h_bound], ...
%        {'Reference Line', ...
%         'Start', ...
%         'End', ...
%         'Bounds'});
% xlabel('x[m]')
% ylabel('y[m]')
% zlabel('z[m]')