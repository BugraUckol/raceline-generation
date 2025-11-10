% F = k1 * (w1 + w2 + w3 + w4)
% Mx = k1 * lp * (w1 + w2 - w3 -w4) 
% My = k1 * lp * (-w1 + w2 + w3 - w4) 
% Mz = k2 * (w1 - w2 + w3 -w4) 

syms k1 k2 lp

w2FM = [k1, k1, k1, k1;
        k1*lp, k1*lp, -k1*lp, -k1*lp;
        -k1*lp, k1*lp, k1*lp, -k1*lp;
        k2, -k2, k2, -k2];

FM2w = w2FM^-1;

% w 0-1
% 5N = 4k1 -> k1 = 1.25
% lp = 0.1061

% k1 =
% 
%     1.2500
% 
% 
% k2 =
% 
%     0.0265
% 
% 
% d =
% 
%     0.2122


k1s = 1.25;
lps = 0.2122;
k2s = 0.0530;

FM2w_num = subs(FM2w, k1, k1s);
FM2w_num = subs(FM2w_num, lp, lps);
FM2w_num = double(subs(FM2w_num, k2, k2s));

w2FM_num = FM2w_num^(-1);

T_max = (w2FM_num * [1,1,1,1]')' * [1,0,0,0]'
Mx_max = (w2FM_num * [1,1,0,0]')' * [0,1,0,0]'
My_max = (w2FM_num * [0,1,1,0]')' * [0,0,1,0]'
Mz_max = (w2FM_num * [1,0,1,0]')' * [0,0,0,1]'

% w2FM_num = subs(w2FM, k1, 0.25);
% w2FM_num = subs(w2FM_num, lp, 10*0.1061);
% w2FM_num = double(subs(w2FM_num, k2, 1*0.1061));








