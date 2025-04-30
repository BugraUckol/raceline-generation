%% Load mat file
% clear, close all, clc
load yas_marina_xy.mat
% load results_reparameterized
num_of_points = 8000;
%% Pre-spline processing
% Close the loop
YasMarina = [YasMarina; YasMarina(1,:)];

% Path calculations
x_y_yaw_curv_s = points_to_yaw_curv_s(YasMarina, 3);

% Shift the track limits
left_side = lateral_shift(x_y_yaw_curv_s,YasMarina(:,3));
right_side = lateral_shift(x_y_yaw_curv_s,-YasMarina(:,4));

% Create final object
x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr = [x_y_yaw_curv_s, left_side(:,1), left_side(:,2), ...
    right_side(:,1), right_side(:,2), YasMarina(:,3), -YasMarina(:,4)];

%% Fit splines
% Fit a spline to center line
spline_x = spline(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,1));
spline_y = spline(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,2));
% Fit splines to road limits
spline_left = spline(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,10));
spline_right = spline(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,11));

t_fine = linspace(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(1,5), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(end,5), num_of_points);
x_smooth = ppval(spline_x, t_fine);
y_smooth = ppval(spline_y, t_fine);
left_smooth = ppval(spline_left, t_fine);
right_smooth = ppval(spline_right, t_fine);

YasMarinaSpline = [x_smooth', y_smooth', left_smooth', right_smooth'];

%% Post-spline processing

% Path calculations
x_y_yaw_curv_s = points_to_yaw_curv_s(YasMarinaSpline, 5);

% Recalculate yaw
% Arclength parameterization with new splines
ppx_s = spline(x_y_yaw_curv_s(:,5), x_smooth);
ppy_s = spline(x_y_yaw_curv_s(:,5), y_smooth);
% p = 0.9;
% [ppx_s, ~] = csaps(x_y_yaw_curv_s(:,5), x_smooth, p);
% [ppy_s, ~] = csaps(x_y_yaw_curv_s(:,5), y_smooth, p);

% Derivatives
dx_ds = ppval(fnder(ppx_s,1), x_y_yaw_curv_s(:,5));
dy_ds = ppval(fnder(ppy_s,1), x_y_yaw_curv_s(:,5));
d2x_ds2 = ppval(fnder(ppx_s,2), x_y_yaw_curv_s(:,5));
d2y_ds2 = ppval(fnder(ppy_s,2), x_y_yaw_curv_s(:,5));
% Yaw
spline_yaw = atan2(dy_ds, dx_ds);
% x_y_yaw_curv_s(:,3) = spline_yaw';

% Recalculate curvature
spline_curvature = (dx_ds .* d2y_ds2 - dy_ds .* d2x_ds2) ./ (dx_ds.^2 + dy_ds.^2).^(3/2);
% x_y_yaw_curv_s(:,4) = spline_curvature';

% Shift the track limits
spline_left_side = lateral_shift(x_y_yaw_curv_s,YasMarinaSpline(:,3));
spline_right_side = lateral_shift(x_y_yaw_curv_s,YasMarinaSpline(:,4));

% Create final object
spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr = [x_y_yaw_curv_s, ... 
    spline_left_side(:,1), spline_left_side(:,2), ...
    spline_right_side(:,1), spline_right_side(:,2), ...
    YasMarinaSpline(:,3), YasMarinaSpline(:,4)];

%% Save
save yas_marina_8000 spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr

% Plot
figure(1)
subplot(1,2,1)
hold on
plot(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,4))
title('Curvature')
grid minor
% for i = 1:max(size(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5)))
%     plot(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(i,5), ...
%         spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(i,4), ...
%         'o', 'DisplayName', num2str(i), 'MarkerFaceColor', 'k', ...
%         'MarkerEdgeColor', 'k', 'MarkerSize', 3);
% end
subplot(1,2,2)
hold on
plot(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,1), spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,2), 'LineWidth', 2);
plot(spline_left_side(:,1), spline_left_side(:,2), 'LineWidth', 2);
plot(spline_right_side(:,1), spline_right_side(:,2), 'LineWidth', 2);
% for i = 1:max(size(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5)))
%     plot(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(i,1), ...
%         spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(i,2), ...
%         'o', 'DisplayName', num2str(i), 'MarkerFaceColor', 'k', ...
%         'MarkerEdgeColor', 'k', 'MarkerSize', 3);
% end
grid minor
title('Map')
daspect([1,1,1])
xlim([min(min(right_side(:,1)), min(right_side(:,2))),...
    max(max(right_side(:,1)), max(right_side(:,2)))])
ylim([min(min(right_side(:,1)), min(right_side(:,2))),...
    max(max(right_side(:,1)), max(right_side(:,2)))])

figure(2)
hold on
plot(x_smooth, y_smooth);
hold on;
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,1), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,2));
hold on
daspect([1,1,1])

figure(3)
subplot(2,1,1)
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,3))
hold on
plot(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), ...
    spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,3))
subplot(2,1,2)
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), ...
    x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,4))
hold on
plot(spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,5), ...
    spline_x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,4))