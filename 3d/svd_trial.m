a = sym("a","real");
b = sym("b","real");
c = sym("c","real");

t = [a b c];

[vec, val] = eig(t' * t)

tr = [0 1 0]';

[vec, val] = eig(tr * tr')

[u,s,v] = svd(tr)