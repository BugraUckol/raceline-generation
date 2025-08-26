clc, clear, close all
% Parameters
nPoints   = 500;          % number of sample points
x         = linspace(-2, 2, nPoints);
y         = 0.6*x.*x;

% First and second derivatives (for curvature and tangent)
dx  = gradient(x);
dy  = gradient(y, x);    % dy/dx
ddx = gradient(dx, x);
ddy = gradient(dy, x);

% Tangent vector (not normalized yet)
Tx = ones(size(x));
Ty = dy;

% Normalize tangent vector
Tnorm = sqrt(Tx.^2 + Ty.^2);
Tx = Tx ./ Tnorm;
Ty = Ty ./ Tnorm;

% Normal vector (perpendicular to tangent)
Nx = -Ty;
Ny = Tx;

% Curvature (kappa = |x'y'' - y'x''| / ( (x'^2+y'^2)^(3/2) ))
xprime  = gradient(x);
yprime  = gradient(y, x);
xpp     = gradient(xprime, x);
ypp     = gradient(yprime, x);
kappa   = abs(xprime.*ypp - yprime.*xpp) ./ ((xprime.^2 + yprime.^2).^(3/2));

% Plotting
figure; hold on; axis equal
plot(x, y, 'k', 'LineWidth', 2, 'LineStyle', '--');  % sine curve

% Parameters for deformed grid
nOffsets = 15;       % number of equidistant offset curves on each side
dOffset  = 0.1;     % spacing between curves
nNormals = 40;      % number of normal lines to draw

% Draw equidistant curves (parallel curves along normals)
for j = -nOffsets:nOffsets
    if j ~= 0
        Xoffset = x + j*dOffset*Nx;
        Yoffset = y + j*dOffset*Ny;
        plot(Xoffset, Yoffset, 'k', 'LineWidth', 0.5)
    end
end

% Draw normal lines at selected points
idx = round(linspace(1, nPoints, nNormals));
for k = idx
    % pick normal line length = max offset
    Xline = x(k) + (-nOffsets*dOffset:nOffsets*dOffset)*Nx(k);
    Yline = y(k) + (-nOffsets*dOffset:nOffsets*dOffset)*Ny(k);
    plot(Xline, Yline, 'k', 'LineWidth', 0.5)
end

% Plot curvature as color on curve (optional)
title('Curvilinear Grid around a Sine Curve')
xlabel('x'); ylabel('y')
