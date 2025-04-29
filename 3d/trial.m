% Infinity-shaped 3D space curve (non-self-intersecting)
% Author: OpenAI / ChatGPT
% Date: 2025-04-30

clear; clc; close all;

% Parameter
t = linspace(0, 2*pi, 1000);  % parametric variable

% Parameters to control shape
a = 1;      % major amplitude (horizontal)
b = 0.5;    % vertical amplitude
c = 0.3;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = a * sin(t);
y = b * sin(2*t);
z = c * cos(t);

% Plot
figure;
plot3(x, y, z, 'k', 'LineWidth', 2);
grid on; axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('3D Infinity Curve (Non-Intersecting)');
view(135, 30); % adjust view angle for clarity