syms alpha beta alpha_dot beta_dot p q r kappa tau s_dot eephi eethe eepsi T v m

cb2w = CB2W(alpha, beta);
cp2w = CE2B([eephi, eethe, eepsi]);

w_w2e_b = [p - beta_dot * sin(alpha);
       q - alpha_dot;
       r + beta_dot * cos(alpha)];

w_w2e_w = cb2w * w_w2e_b;

w_p2e_p = s_dot * [tau; 0; kappa];
w_p2e_w = cp2w * w_p2e_p;

w_w2p_w = w_w2e_w - w_p2e_w;

% beta alpha dots
adot = ((q*cos(beta) - p*cos(alpha)*sin(beta) - r*sin(alpha)*sin(beta)) + T*cos(alpha)/m/v) / cos(beta);
bdot = (-r*cos(alpha) + p*sin(alpha) - sin(alpha)*sin(beta)*T/m/v);

% subbed

w_w2p_w = simplify(subs(w_w2p_w, alpha_dot, adot));
w_w2p_w = simplify(subs(w_w2p_w, beta_dot, bdot));

% sdot
syms sdot v een
sdot = (v * cos(eepsi) * cos(eethe)) / (1 - kappa * een);

% subbed
w_w2p_w = simplify(subs(w_w2p_w, s_dot, sdot));

E = [eephi, eethe, eepsi];
ed_w2p_w = simplify(W2ED(w_w2p_w, E));

% spatial derivative
w_w2p_w = simplify(ed_w2p_w / sdot)














