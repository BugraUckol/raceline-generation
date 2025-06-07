function [ED] = W2ED(W,E)

phi = E(1);
the = E(2);

H = [1, sin(phi) * tan(the), cos(phi) * tan(the);
    0,             cos(phi),           -sin(phi);
    0,  sin(phi) * sec(the), cos(phi) * sec(the)];

ED = H * W;

end