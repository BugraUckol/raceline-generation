function [FM2w_num] = createAllocationMatrix(T_max, Mxy_max, Mz_max)

% F = k1 * (w1 + w2 + w3 + w4)
% Mx = k1 * lp * (w1 + w2 - w3 -w4) 
% My = k1 * lp * (-w1 + w2 + w3 - w4) 
% Mz = k2 * (w1 - w2 + w3 -w4) 

syms k1 k2 lp

w2FM = [k1, k1, k1, k1;
        k1*lp, k1*lp, -k1*lp, -k1*lp;
        -k1*lp, k1*lp, k1*lp, -k1*lp;
        k2, -k2, k2, -k2];

k1_num = T_max / 4
k2_num = Mz_max / 4
lp_num = Mxy_max / 4 / k1_num;
l = lp_num * 2 / sqrt(2)

w2FM_num = double(subs(subs(subs(w2FM,k1,k1_num), k2,k2_num), lp, lp_num));

FM2w_num = w2FM_num^-1;

end






