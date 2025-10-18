function dcm_new = dcmUpdateRk4(dcm,omega,ts)

% Function Definitions
f1 = @(x1, w) skew(w) * x1;

% 4th Order Explicit Runge Kutta Integration
dcm_new = RK4(dcm, omega, f1, ts);

end
