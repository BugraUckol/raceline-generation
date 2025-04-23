% Load mat file
clc, clear, close all
load yas_marina_xy.mat

% Close the loop
YasMarina = [YasMarina; YasMarina(1,:)];

% Path calculations
x_y_yaw_curv_s = points_to_yaw_curv_s(YasMarina, 1);

% Shift the track limits
left_side = lateral_shift(x_y_yaw_curv_s,YasMarina(:,3));
right_side = lateral_shift(x_y_yaw_curv_s,-YasMarina(:,4));

% Create final object
x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr = [x_y_yaw_curv_s, left_side(:,1), left_side(:,2), ...
    right_side(:,1), right_side(:,2), YasMarina(:,3), -YasMarina(:,4)];

% Save
save yas_marina_full_narrow_window x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr

% Plot
figure(1)
subplot(1,2,1)
hold on
plot(x_y_yaw_curv_s(:,5), x_y_yaw_curv_s(:,4))
title('Curvature')
grid minor
subplot(1,2,2)
hold on
plot(x_y_yaw_curv_s(:,1), x_y_yaw_curv_s(:,2), 'LineWidth', 2);
plot(left_side(:,1), left_side(:,2), 'LineWidth', 2);
plot(right_side(:,1), right_side(:,2), 'LineWidth', 2);
grid minor
title('Map')
daspect([1,1,1])
xlim([min(min(right_side(:,1)), min(right_side(:,2))),...
    max(max(right_side(:,1)), max(right_side(:,2)))])
ylim([min(min(right_side(:,1)), min(right_side(:,2))),...
    max(max(right_side(:,1)), max(right_side(:,2)))])

figure(2)
hold on
plot(x_y_yaw_curv_s_lx_ly_rx_ry_dl_dr(:,1), x_y_yaw_curv_s(:,2), 'LineWidth', 2);
hold on
