%% Prep
clc, clear, close all

%% Symbols Calculations
m = sym("m","real");
I_xx = sym("I_xx","real");
I_yy = sym("I_yy","real");
I_zz = sym("I_zz","real");
Cd = sym("Cd","real");

t = sym("t","real");

u = sym("u","real");
v = sym("v","real");
w = sym("w","real");

p = sym("p","real");
q = sym("q","real");
r = sym("r","real");

s = sym("s","real");
en = sym("en","real");
eb = sym("eb","real");

u_dot = sym("u_dot","real");
v_dot = sym("v_dot","real");
w_dot = sym("w_dot","real");

p_dot = sym("p_dot","real");
q_dot = sym("q_dot","real");
r_dot = sym("r_dot","real");

s_dot = sym("s_dot","real");
en_dot = sym("en_dot","real");
eb_dot = sym("eb_dot","real");

e_phi = sym("e_phi","real");
e_the = sym("e_the","real");
e_psi = sym("e_psi","real");

e_phi_dot = sym("e_phi_dot","real");
e_the_dot = sym("e_the_dot","real");
e_psi_dot = sym("e_psi_dot","real");

tau = sym("tau","real");
kappa1 = sym("kappa1","real");
kappa2 = sym("kappa2","real");

Mx = sym("Mx","real");
My = sym("My","real");
Mz = sym("Mz","real");
Fz = sym("Fz","real");

gx_p = sym("gx_p","real");
gy_p = sym("gy_p","real");
gz_p = sym("gz_p","real");

%% Definitions
% Transformations
E = [e_phi, e_the, e_psi];
Cb2p = CB2E(E);
Cp2b = CE2B(E);

% Vectors
w_be_b = [p; q; r];
w_be_p = Cb2p * w_be_b;
w_pe_p = s_dot * [tau; kappa1; kappa2];
w_pe_b = Cp2b * w_pe_p;
R_bp_p = [0; en; eb];
De_R_pe_p = [s_dot; 0; 0];
Dp_R_bp_p = [0; en_dot; eb_dot];
De_R_be_b = [u; v; w];
De_R_be_p = Cb2p * De_R_be_b;
Db_w_be_b = [p_dot; q_dot; r_dot];
E_dot = [e_phi_dot; e_the_dot; e_psi_dot];
F_t = [0; 0; Fz];
F_d = -Cd * De_R_be_b;
F_b = F_t + F_d;
M_b = [Mx; My; Mz];
I_b = [I_xx, 0, 0; 0, I_yy, 0; 0, 0, I_zz];
V_dot = [u_dot; v_dot; w_dot];
W_dot = [p_dot; q_dot; r_dot];
g_p = [gx_p; gy_p; gz_p];

%% Newton-Euler Equations
eq_newton = V_dot == F_b/m + Cp2b * g_p - cross(w_be_b, De_R_be_b);
eq_euler = W_dot == I_b^-1 * (M_b - cross(w_be_b, I_b * w_be_b));

%% Kinematical Loop Closure
% Position of the UAV wrt. Earth = Position of the UAV wrt. Path + Position of the Path wrt. Earth
% R_be = R_bp + R_pe

% Differentiating the expression above in Earth
% De_R_be = De_R_bp + De_R_pe

% Using Corilois Transport Theorem on The De_R_pe. Also De_R_bp = V
% V = De_R_pe + Dp_R_bp + w_pe x R_bp

% Resolving the equation in Fp
eq_loop = Cb2p * De_R_be_b == De_R_pe_p + Dp_R_bp_p + cross(w_pe_p, R_bp_p);

% Isolating state derivatives
eq_s_dot = rhs((isolate(eq_loop(1), s_dot)));
eq_en_dot = rhs((isolate(eq_loop(2), en_dot)));
eq_eb_dot = rhs((isolate(eq_loop(3), eb_dot)));

%% Attitude Dynamics
% w_bp_p == w_be_p - w_pe_p;
eq_attitude_dot = W2ED(w_be_b - w_pe_b, E);

%% Spatial Equations
eq_t_prime = eq_s_dot^-1;
eq_en_prime = simplify(subs(eq_en_dot, s_dot, eq_s_dot) * eq_t_prime);
eq_eb_prime = simplify(subs(eq_eb_dot, s_dot, eq_s_dot) * eq_t_prime);

eq_attitude_spatial = simplify(subs(eq_attitude_dot, s_dot, eq_s_dot) * eq_t_prime);

eq_newton_spatial = simplify(rhs(eq_newton) * eq_t_prime);
eq_euler_spatial = simplify(rhs(eq_euler) * eq_t_prime);

[eq_t_prime; eq_en_prime; eq_eb_prime; eq_attitude_spatial;
    eq_newton_spatial; eq_euler_spatial]

































