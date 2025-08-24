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

X = [t, ey, ez, e_psi, e_the, e_phi, v, alpha, beta]
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
pqr_lim = 20;
T_max = 5;
T_min = 0;
m = 0.250;
cv_lin = 0.00735;

%% Initial conditions
t0 = 0;
e_n0 = 0.0;
e_b0 = 0.0;
e_phi0 = 0.0;
e_psi0 = 0.0;
e_the0 = 0.0;
v0 = 10.0;
alpha0 = 0.0;
beta0 = 0.0;

%% Setting Optimization Problem
size_vec = floor(size(path.s_arr));
N = size_vec(2) - 1;

opti = casadi.Opti();

X = opti.variable(9, N+1); % state trajectory in path frame
t = X(1,:);
en = X(2,:);
eb = X(3,:);
ephi = X(4,:);
ethe = X(5,:);
epsi = X(6,:);
vv = X(7,:);
alphav = X(8,:);
betav = X(9,:);

U = opti.variable(4,N);   % Angular rates
T_com = U(1,:);
p_com = U(2,:);
q_com = U(3,:);
r_com = U(4,:);

% Cost Function
% Minimize angular velocity inputs p and q
% opti.minimize(1.0 * U(1,:) * U(1,:)' + 1.0 * U(2,:) * U(2,:)');
% opti.minimize(1.0 * (en * en') + 1.0 * (eb * eb'));
% opti.minimize(1.0 * (t * t'));
opti.minimize(t(end));


% System dynamics
% x* = [t, ey, eb, e_phi, e_the, e_psi, v, alpha, beta] states in spatial formulation
% u* = [T, p, q, r] inputs in spatial formulation
% a* = [(-tan(beta) * (q*cos(beta) - p*cos(alpha)*sin(beta) - r*sin(alpha)*sin(beta) + T * cos(alpha) / m / v) + q*sin(beta) + p*cos(alpha)*cos(beta) + r*cos(beta)*sin(alpha)),(T * cos(alpha) / m / v),(-T * sin(alpha) * sin(beta) / m / v),kappa,tau] auxilary inputs for better calculation
f = @(tt,een,eeb,eephi,eethe,eepsi,v,alpha,beta,T,p,q,r,kappa,tau) [
    % dt / ds
    (1 - kappa * een) / (v * cos(eepsi) * cos(eethe));
    % dey / ds
    (1 - kappa * een) * tan(eepsi) + tau * eeb;
    % deb / ds
    -(1 - kappa * een) * tan(eethe) / cos(eepsi) - tau * een;
    % dephi / ds (phi of the of the velocity frame wrt. the path frame)
    -(T*cos(alpha)*cos(eethe)*sin(beta) - m*p*v*cos(alpha)*cos(eethe) - m*r*v*cos(eethe)*sin(alpha) + T*cos(alpha)*cos(beta)*sin(eephi)*sin(eethe) + T*cos(beta)*cos(eephi)*sin(alpha)*sin(beta)*sin(eethe) + m*tau*v^2*cos(beta)*cos(eepsi)^2*cos(eethe) - T*een*kappa*cos(alpha)*cos(eethe)*sin(beta) + een*kappa*m*p*v*cos(alpha)*cos(eethe) + een*kappa*m*r*v*cos(eethe)*sin(alpha) - T*een*kappa*cos(alpha)*cos(beta)*sin(eephi)*sin(eethe) - T*een*kappa*cos(beta)*cos(eephi)*sin(alpha)*sin(beta)*sin(eethe))/(m*v^2*cos(beta)*cos(eepsi)*cos(eethe)^2)
    % dethe / ds (theta of the of the velocity frame wrt. the path frame)
    (m*tau*cos(eepsi)*cos(eethe)*sin(eepsi)*v^2 - T*cos(alpha)*cos(eephi) + T*sin(alpha)*sin(beta)*sin(eephi) + T*een*kappa*cos(alpha)*cos(eephi) - T*een*kappa*sin(alpha)*sin(beta)*sin(eephi))/(m*v^2*cos(eepsi)*cos(eethe))
    % depsi / ds (psi of the of the velocity frame wrt. the path frame)
    -(T*cos(alpha)*sin(eephi) + T*cos(eephi)*sin(alpha)*sin(beta) - T*een*kappa*cos(alpha)*sin(eephi) + kappa*m*v^2*cos(eepsi)*cos(eethe)^2 + m*tau*v^2*cos(eepsi)^2*cos(eethe)*sin(eethe) - T*een*kappa*cos(eephi)*sin(alpha)*sin(beta))/(m*v^2*cos(eepsi)*cos(eethe)^2)
    % dv / ds
    (T*sin(alpha)*cos(beta) / m) * ((1 - kappa * een) / (v * cos(eepsi) * cos(eethe)));
    % dalpha / ds
    ((q*cos(beta) - p*cos(alpha)*sin(beta) - r*sin(alpha)*sin(beta)) + T*cos(alpha)/m/v) / cos(beta) * ((1 - kappa * een) / (v * cos(eepsi) * cos(eethe)));
    % dbeta / ds
    (-r*cos(alpha) + p*sin(alpha) - sin(alpha)*sin(beta)*T/m/v) * ((1 - kappa * een) / (v * cos(eepsi) * cos(eethe)));
    ];

for k=1:N % loop over control intervals
   kappa = path.kappa_arr(k);
   tau = path.tau_arr(k);
   ds = path.s_arr(k + 1) ...
       - path.s_arr(k);

   % 1st Order Explicit Euler's Integration
   k1 = f(X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), X(6,k),...
       X(7,k), X(8,k), X(9,k), U(1,k), U(2,k), U(3,k), U(4,k), ...
       kappa, tau);
   x_next = X(:,k) + ds * k1;
   
   % Multiple shooting
   opti.subject_to(X(:,k+1)==x_next); % close the gaps

end
opti.subject_to(t(2:N+1) > t(1:N)) % Time must increase!

% Constraints
opti.subject_to(U(1,:) <= T_max);
opti.subject_to(U(2,:) <= pqr_lim);
opti.subject_to(U(3,:) <= pqr_lim);
opti.subject_to(U(4,:) <= pqr_lim);
opti.subject_to(U(1,:) >= T_min);
opti.subject_to(U(2,:) >= -pqr_lim);
opti.subject_to(U(3,:) >= -pqr_lim);
opti.subject_to(U(4,:) >= -pqr_lim);

% opti.subject_to(X(4,:) <= 0.2);
% opti.subject_to(X(4,:) >= 0.2);
% opti.subject_to(X(5,:) <= 0.5);
% opti.subject_to(X(5,:) >= -0.5);
% opti.subject_to(X(6,:) <= 0.5);
% opti.subject_to(X(6,:) >= -0.5);
% 
% opti.subject_to(U(1,:) == 0.0);
% opti.subject_to(X(4,:) >= -0.1);

opti.subject_to(X(2,:) <= 0.5);
opti.subject_to(X(2,:) >= -0.5);
opti.subject_to(X(3,:) <= 0.5);
opti.subject_to(X(3,:) >= -0.5);
% opti.subject_to(X(3,:) == 0.0);

opti.subject_to(X(4,:) <= pi);
opti.subject_to(X(4,:) >= -pi);
opti.subject_to(X(5,:) <= pi);
opti.subject_to(X(5,:) >= -pi);
opti.subject_to(X(6,:) <= pi);
opti.subject_to(X(6,:) >= -pi);

opti.subject_to(X(8,:) <= pi);
opti.subject_to(X(8,:) >= -pi);
opti.subject_to(X(9,:) <= pi);
opti.subject_to(X(9,:) >= -pi);

% opti.subject_to(X(7,:) >= 7.8);

opti.subject_to(t(1) == 0.0);
% opti.subject_to(en(1) == e_n0);
% opti.subject_to(eb(1) == e_b0);
% opti.subject_to(ephi(1) == e_phi0);
% opti.subject_to(ethe(1) == e_the0);
% opti.subject_to(epsi(1) == e_psi0);
% opti.subject_to(vv(1) == v0);
% opti.subject_to(alphav(1) == alpha0);
% opti.subject_to(betav(1) == beta0);
opti.subject_to(en(1) == en(end));
opti.subject_to(eb(1) == eb(end));
opti.subject_to(ephi(1) == ephi(end));
opti.subject_to(ethe(1) == ethe(end));
opti.subject_to(epsi(1) == epsi(end));
opti.subject_to(vv(1) == vv(end));
opti.subject_to(alphav(1) == alphav(end));
opti.subject_to(betav(1) == betav(end));

% prevsol = load('prevsol_circle');
% opti.set_initial(t, prevsol.time_arr);
% opti.set_initial(en, prevsol.en_arr);
% opti.set_initial(eb, prevsol.eb_arr);
% opti.set_initial(ephi, prevsol.ephi_arr);
% opti.set_initial(ethe, prevsol.ethe_arr);
% opti.set_initial(epsi, prevsol.epsi_arr);
% opti.set_initial(vv, prevsol.v_arr);
% opti.set_initial(alphav, prevsol.alpha_arr);
% opti.set_initial(betav, prevsol.beta_arr);
% 
% opti.set_initial(T_com, prevsol.T_com_arr);
% opti.set_initial(p_com, prevsol.p_com_arr);
% opti.set_initial(q_com, prevsol.q_com_arr);
% opti.set_initial(r_com, prevsol.r_com_arr);
opti.set_initial(vv, 10);

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
v_arr = sol.value(vv);
alpha_arr = sol.value(alphav);
beta_arr = sol.value(betav);

T_com_arr = sol.value(T_com);
p_com_arr = sol.value(p_com);
q_com_arr = sol.value(q_com);
r_com_arr = sol.value(r_com);

save prevsol_circle time_arr en_arr eb_arr ephi_arr ethe_arr epsi_arr v_arr alpha_arr beta_arr T_com_arr p_com_arr q_com_arr r_com_arr

%% 3D Recreation
x_arr = zeros(1,length(distance_arr));
y_arr = zeros(1,length(distance_arr));
z_arr = zeros(1,length(distance_arr));

k = 0;

edges = [0, 0, 0, 0, 0;
         0.5, -0.5, -0.5, 0.5, 0.5;
         0.5, 0.5, -0.5, -0.5, 0.5];
edges_e_list = [];
p0 = [path.x_arr(1); path.y_arr(1); path.z_arr(1)];
dcm_v_p = CB2E([e_phi0, e_the0, e_psi0]);
dcm_p_e = CB2E([path.roll_arr(1), path.pitch_arr(1), path.yaw_arr(1)]);
dcm_b_v = CB2W(alpha_arr(1), beta_arr(1));
r = dcm_p_e * [0; e_n0; e_b0] + p0;
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

    dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
    dcm_b_e = CB2E([ephi_arr(k), ephi_arr(k), ephi_arr(k)]);
    pt = [path.x_arr(k); path.y_arr(k); path.z_arr(k)];
    ptx = [pt, pt + dcm_p_e(:,1) * 1.5];
    pty = [pt, pt + dcm_p_e(:,2) * 1.5];
    ptz = [pt, pt + dcm_p_e(:,3) * 1.5];
    % plot3(ptx(1,:), ptx(2,:), ptx(3,:), 'LineWidth', 3, 'Color', 'r');
    % plot3(pty(1,:), pty(2,:), pty(3,:), 'LineWidth', 3, 'Color', 'g');
    % plot3(ptz(1,:), ptz(2,:), ptz(3,:), 'LineWidth', 3, 'Color', 'b');

    edges_e = pt + dcm_p_e * edges;
    edges_e_list = cat(3, edges_e_list, edges_e);
end

%% Plots
setDefaultFigureProperties()
set(0,'DefaultFigureWindowStyle','docked')
figure(1)
subplot(3,2,1)
plot(distance_arr, time_arr);
grid minor
title('Curvilinear Distance vs Time')
xlabel('Curvilinear Distance [m]')
hold on

subplot(3,2,2)
plot(distance_arr, en_arr), hold on
plot(distance_arr, eb_arr)
title('Curvilinear Distance vs Cartesian Deviations')
legend('e_n', 'e_b')
xlabel('Curvilinear Distance [m]')
ylabel('Cartesian Deviation [m]')
grid minor

subplot(3,2,3)
plot(distance_arr, ephi_arr), hold on
plot(distance_arr, ethe_arr)
plot(distance_arr, epsi_arr)
title(['Curvilinear Distance vs Euler Angles of $$\mathcal{F_V}$$ wrt. $$\mathcal{F_P}$$']),...
    legend('\gamma_{ex}', '\gamma_{ey}', '\gamma_{ez}')
xlabel('Curvilinear Distance [m]')
ylabel('Angle [rad]')
grid minor

subplot(3,2,4)
plot(distance_arr, v_arr)
hold on
title('Curvilinear Distance vs Velocity')
xlabel('Curvilinear Distance [m]')
ylabel('Velocity [m/s]')
grid minor

subplot(3,2,5)
plot(distance_arr, alpha_arr)
hold on
title('Curvilinear Distance vs Angle of Attack ($$\alpha$$)')
xlabel('Curvilinear Distance [m]')
ylabel('Angle [rad]')
grid minor

subplot(3,2,6)
plot(distance_arr, beta_arr)
hold on
title('Curvilinear Distance vs Angle of Sideslip ($$\beta$$)')
xlabel('Curvilinear Distance [m]')
ylabel('Angle [rad]')
grid minor

%%
figure(2)
subplot(2,1,1)
plot(distance_arr(1:end-1), T_com_arr); hold on;
title('Curvilinear Distance vs Thrust')
xlabel('Curvilinear Distance [m]')
ylabel('Thrust [N]')
grid minor
subplot(2,1,2)
plot(distance_arr(1:end-1), p_com_arr); hold on;
plot(distance_arr(1:end-1), q_com_arr);
plot(distance_arr(1:end-1), r_com_arr);
grid minor
title('Curvilinear Distance vs Angular Rates')
xlabel('Curvilinear Distance [m]')
ylabel('Angular Velocity [Rad/s]')
legend('p','q','r')
%%
figure(3)
plot3(path.x_arr, path.y_arr, path.z_arr, 'LineWidth', 2); hold on
plot3(x_arr, y_arr, z_arr, 'LineWidth', 2);
daspect([1,1,1])
xlabel('East'), ylabel('North'), zlabel('Up'), title('Trefoil Tube and The Optimal Solution')
grid minor; box on
p0 = [path.x_arr(1); path.y_arr(1); path.z_arr(1)];
dcm_v_p = CB2E([e_phi0, e_the0, e_psi0]);
dcm_p_e = CB2E([path.roll_arr(1), path.pitch_arr(1), path.yaw_arr(1)]);
dcm_b_v = CB2W(alpha_arr(1), beta_arr(1));
r = dcm_p_e * [0; e_n0; e_b0] + p0;
% dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
dcm_v_e  = dcm_p_e * dcm_v_p;
r_hist = [r];
for k=1:1:N-1
    cla;
    dcm_p_e = CB2E([path.roll_arr(k), path.pitch_arr(k), path.yaw_arr(k)]);
    plot3(path.x_arr, path.y_arr, path.z_arr, 'LineWidth', 2)
    plot3(x_arr, y_arr, z_arr, 'LineWidth', 2)
    % plot3(r_hist(1,:), r_hist(2,:), r_hist(3,:), 'LineWidth', 3)
    
    px = [r, r + dcm_v_e(:,1)];
    py = [r, r + dcm_v_e(:,2)];
    pz = [r, r + dcm_v_e(:,3)];
    
    % plot3(px(1,:), px(2,:), px(3,:), 'LineWidth', 2, 'Color', 'r');
    % plot3(py(1,:), py(2,:), py(3,:), 'LineWidth', 2, 'Color', 'g');
    % plot3(pz(1,:), pz(2,:), pz(3,:), 'LineWidth', 2, 'Color', 'b');

    dt = time_arr(k+1) - time_arr(k);
    ds = distance_arr(k+1) - distance_arr(k);

    r = r + dcm_v_e * [v0; 0; 0] * dt;

    r_hist = [r_hist, r];

    w_be = [p_com_arr(k), q_com_arr(k), r_com_arr(k)];

    eilist_size = size(edges_e_list);
    for ei = 2:20:eilist_size(3)-1
    edges_ei = edges_e_list(:,:,ei);
    plot3(edges_ei(1,:), edges_ei(2,:), edges_ei(3,:), 'k', 'LineWidth', 0.15), hold on
    end
    edges_ei = edges_e_list(:,:,1);
    plot3(edges_ei(1,:), edges_ei(2,:), edges_ei(3,:), 'g', 'LineWidth', 4), hold on
    edges_ei = edges_e_list(:,:,end-1);
    plot3(edges_ei(1,:), edges_ei(2,:), edges_ei(3,:), 'r', 'LineWidth', 4), hold on

    plot3(squeeze(edges_e_list(1,1,:)), squeeze(edges_e_list(2,1,:)), squeeze(edges_e_list(3,1,:)), 'k', 'LineWidth', 1.0)
    plot3(squeeze(edges_e_list(1,2,:)), squeeze(edges_e_list(2,2,:)), squeeze(edges_e_list(3,2,:)), 'k', 'LineWidth', 1.0)
    plot3(squeeze(edges_e_list(1,3,:)), squeeze(edges_e_list(2,3,:)), squeeze(edges_e_list(3,3,:)), 'k', 'LineWidth', 1.0)
    plot3(squeeze(edges_e_list(1,4,:)), squeeze(edges_e_list(2,4,:)), squeeze(edges_e_list(3,4,:)), 'k', 'LineWidth', 1.0)
    % plot3(edges_e(1,:), edges_e(2,:), edges_e(3,:), 'k', 'LineWidth', 0.5)
     

    % [spx, spy, spz] = sphere(10);
    % surf(path.x_arr(1) + 0.5*spx, path.y_arr(1) + 0.5*spy, path.z_arr(1) + 0.5*spz, 'FaceAlpha',0.5, 'FaceColor', 'r');

    % Ground truth of velocity frame
    dcm_gt_p = CB2E([ephi_arr(k), ethe_arr(k), epsi_arr(k)]);
    dcm_gt_e = dcm_p_e * dcm_gt_p;
    rgt = [x_arr(k); y_arr(k); z_arr(k)];
    pgtx = [rgt, rgt + dcm_gt_e(:,1) * 1.5];
    pgty = [rgt, rgt + dcm_gt_e(:,2) * 1.5];
    pgtz = [rgt, rgt + dcm_gt_e(:,3) * 1.5];
    
    % Plot Velocity Frame
    plot3(pgtx(1,:), pgtx(2,:), pgtx(3,:), 'LineWidth', 3, 'Color', 'm');
    % plot3(pgty(1,:), pgty(2,:), pgty(3,:), 'LineWidth', 3, 'Color', 'g');
    % plot3(pgtz(1,:), pgtz(2,:), pgtz(3,:), 'LineWidth', 3, 'Color', 'b');

    % Plot Thrust Vector
    dcm_b_v = CB2W(alpha_arr(k), beta_arr(k));
    dcm_b_e = dcm_gt_e * dcm_b_v;
    pvx = [rgt, rgt + dcm_b_e(:,1) * 1];
    pvy = [rgt, rgt + dcm_b_e(:,2) * 1];
    pvz = [rgt, rgt + dcm_b_e(:,3) * 1];
    plot3(pvx(1,:), pvx(2,:), pvx(3,:), 'LineWidth', 3, 'Color', 'r');
    plot3(pvy(1,:), pvy(2,:), pvy(3,:), 'LineWidth', 3, 'Color', 'g');
    plot3(pvz(1,:), pvz(2,:), pvz(3,:), 'LineWidth', 3, 'Color', 'b');

    drawnow;
    daspect([1,1,1])
    pause(0.05)

    if mod(k,10) == 0
        k;
    end
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
