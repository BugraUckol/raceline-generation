%% Prep
clc, clear, close all

%% Declearations
m = sym("m","real");
I_xx = sym("I_xx","real");
I_yy = sym("I_yy","real");
I_zz = sym("I_zz","real");

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
kappa = sym("kappa","real");

Mx = sym("Mx","real");
My = sym("My","real");
Mz = sym("Mz","real");
Fz = sym("Fz","real");

%% Transformations
E = [e_phi, e_the, e_psi];
Cb2p = CB2E(E);
Cp2b = CE2B(E);

%% Definitions
w_be_b = [p; q; r];
w_pe_p = s_dot * [tau; 0; kappa];
R_bp_p = [0; en; eb];
De_R_pe_p = [s_dot; 0; 0];
Dp_R_bp_p = [0; en_dot; eb_dot];
De_R_be_b = [u; v; w];
De_R_be_p = Cb2p * De_R_be_b;
Db_w_be_b = [p_dot; q_dot; r_dot];
E_dot = [e_phi_dot; e_the_dot; e_psi_dot];

%% Rigid Body Dynamics
F_b = [0; 0; Fz];
M_b = [Mx; My; Mz];
I_b = [I_xx, 0, 0; 0, I_yy, 0; 0, 0, I_zz];

V_dot_lhs = [u_dot; v_dot; w_dot];
V_dot_rhs = F_b/m - cross(w_be_b, De_R_be_b);
eq1 = V_dot_lhs == V_dot_rhs; 

attitude_rhs = I_b^-1 * (M_b - cross(w_be_b, I_b * w_be_b));
eq2 = Db_w_be_b == attitude_rhs;
%% Translational Relations
kinematic_lhs = De_R_be_p;
kinematic_rhs = De_R_pe_p + Dp_R_bp_p  + cross(w_pe_p, R_bp_p);
eq3 = kinematic_lhs == kinematic_rhs;

%% Rotational Relations
dsdt = u * cos(e_psi) * cos(e_the) / (1 - kappa * en);
dtds = 1 / dsdt;

w_pe_b = Cp2b * w_pe_p;

w_bp_b = w_be_b - w_pe_b;
w_b2p_b_ds = w_bp_b * dtds;

eq4_rhs = simplify(W2ED(w_bp_b, E));
eq4 = E_dot == eq4_rhs;

%% Equations
eq_u_dot = rhs((isolate(eq1(1), u_dot)));
eq_v_dot = rhs((isolate(eq1(2), v_dot)));
eq_w_dot = rhs((isolate(eq1(3), w_dot)));
eq_p_dot = rhs((isolate(eq2(1), p_dot)));
eq_q_dot = rhs((isolate(eq2(2), q_dot)));
eq_r_dot = rhs((isolate(eq2(3), r_dot)));
eq_s_dot = rhs((isolate(eq3(1), s_dot)));
eq_en_dot = rhs((isolate(eq3(2), en_dot)));
eq_eb_dot = rhs((isolate(eq3(3), eb_dot)));
eq_phi_dot = rhs((isolate(eq4(1), e_phi_dot)));
eq_the_dot = rhs((isolate(eq4(2), e_the_dot)));
eq_psi_dot = rhs((isolate(eq4(3), e_psi_dot)));

%% Spatial transformation
eq_t_ds = 1/eq_s_dot;

eq_u_ds = simplify(eq_u_dot * eq_t_ds);
eq_v_ds = simplify(eq_v_dot * eq_t_ds);
eq_w_ds = simplify(eq_w_dot * eq_t_ds);

eq_p_ds = simplify(eq_p_dot * eq_t_ds);
eq_q_ds = simplify(eq_q_dot * eq_t_ds);
eq_r_ds = simplify(eq_r_dot * eq_t_ds);

eq_en_ds = simplify(subs(eq_en_dot, s_dot, eq_s_dot) * eq_t_ds);
eq_eb_ds = simplify(subs(eq_eb_dot, s_dot, eq_s_dot) * eq_t_ds);

eq_phi_ds = simplify(subs(eq_phi_dot, s_dot, eq_s_dot) * eq_t_ds);
eq_the_ds = simplify(subs(eq_the_dot, s_dot, eq_s_dot) * eq_t_ds);
eq_psi_ds = simplify(subs(eq_psi_dot, s_dot, eq_s_dot) * eq_t_ds);

%% Accelerations
syms t;
syms ut(t);
syms vt(t);
syms wt(t);
syms e_phit(t);
syms e_thet(t);
syms e_psit(t);
syms sdt(t);

eq_s_dot(t) = -(w*(sin(e_phit)*sin(e_psit) + cos(e_phit)*cos(e_psit)*sin(e_thet)) - v*(cos(e_phit)*sin(e_psit) - cos(e_psit)*sin(e_phit)*sin(e_thet)) + u*cos(e_psit)*cos(e_thet))/(en*kappa - 1)

eq_s_dot = subs(eq_s_dot, u, ut);
eq_s_dot = subs(eq_s_dot, v, vt);
eq_s_dot = subs(eq_s_dot, w, wt);

eq_d_ddot = diff(eq_s_dot, t)