f = Function('f', {vertcat(X(:),U(:))}, {opti.f});
g = Function('g', {vertcat(X(:),U(:))}, {vertcat(opti.g)});

objective = @(z) full(f(z));
constraints = @(z) deal([], full(g(z)));

nvars = 9*(N+1) + 4*N;

lb = -inf(nvars,1);
ub = inf(nvars,1);

problem = createOptimProblem('fmincon', 'objective', objective, 'x0', ...
    vertcat(reshape(sol.value(X),99,[]), reshape(sol.value(U),40,[])), ...
    'lb', lb, 'ub', ub, 'nonlcon', constraints, ...
    'options', optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
    'MaxIterations',10000, 'MaxFunctionEvaluations', 500000, 'ConstraintTolerance', 1e-4));

gs = GlobalSearch('NumTrialPoints',500);
[z_global, fval_global] = run(gs, problem);



