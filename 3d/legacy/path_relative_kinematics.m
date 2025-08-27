clc, clear, close all

phi = sym("eephi","real");
the = sym("eethe","real");
psi = sym("eepsi","real");

p = sym("p","real");
q = sym("q","real");
r = sym("r","real");

tau = sym("tau","real");
kappa = sym("kappa","real");

v0 = sym("v0","real");
een = sym("een","real");

E = [phi, the, psi];

Cb2p = CB2E(E);
Cp2b = CE2B(E);

dsdt = v0 * cos(psi) * cos(the) / (1 - kappa * een);
dtds = 1 / dsdt;

w_p2e_p = [tau; 0; kappa] * dsdt;
w_p2e_b = Cp2b * w_p2e_p;

w_b2e_b = [p; q; r];

w_b2p_b = w_b2e_b - w_p2e_b;
w_b2p_b_ds = w_b2p_b * dtds;

E_dot = simplify(W2ED(w_b2p_b, E));
E_ds = simplify(W2ED(w_b2p_b_ds, E))