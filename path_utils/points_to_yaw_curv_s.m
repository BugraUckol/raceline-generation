function [x_y_yaw_curv_s] = points_to_yaw_curv_s(x_y_arr, half_window_size)
%points_to_yaw_and_curv Calculate a path's yaw and curvature
%   This function calculates a given path's yaw and curvature values
%   with numerical methods. A basic window-filtering approach is used.

% Pre allocation
x_y_yaw_curv_s = zeros(length(x_y_arr(:,1)), 5);
x_y_yaw_curv_s_padded = zeros(length(x_y_arr(:,1)) + 2*half_window_size, 2);
arr_length = length(x_y_arr(:,1));

% Filling the padded curve X
x_y_yaw_curv_s_padded(1 : half_window_size, 1) = ...
    x_y_arr(arr_length-half_window_size + 1 : end, 1);
x_y_yaw_curv_s_padded(half_window_size + arr_length + 1 : end, 1) = ...
    x_y_arr(1:half_window_size, 1);
x_y_yaw_curv_s_padded(half_window_size + 1 : arr_length + ...
    half_window_size, 1) = x_y_arr(:,1);

% Filling the oadded curve Y
x_y_yaw_curv_s_padded(1 : half_window_size, 2) = ...
    x_y_arr(arr_length-half_window_size + 1 : end, 2);
x_y_yaw_curv_s_padded(half_window_size + arr_length + 1 : end, 2) = ...
    x_y_arr(1:half_window_size, 2);
x_y_yaw_curv_s_padded(half_window_size + 1 : arr_length + ...
    half_window_size, 2) = x_y_arr(:,2);

% Yaw, curvature and curvilinear distance calculations
for idx = 1:arr_length
    % Points
    padded_idx = idx + half_window_size;
   
    p1x = x_y_yaw_curv_s_padded(padded_idx - half_window_size, 1);
    p1y = x_y_yaw_curv_s_padded(padded_idx - half_window_size, 2);

    p2x = x_y_yaw_curv_s_padded(padded_idx, 1);
    p2y = x_y_yaw_curv_s_padded(padded_idx, 2);

    p3x = x_y_yaw_curv_s_padded(padded_idx + half_window_size, 1);
    p3y = x_y_yaw_curv_s_padded(padded_idx + half_window_size, 2);

    % Yaw calculation
    yaw = atan2((p3y - p1y), (p3x - p1x));

    % Curvature calculation
    dist1 = sqrt((p1x - p2x)^2 + (p1y - p2y)^2);
    dist2 = sqrt((p1x - p3x)^2 + (p1y - p3y)^2);
    dist3 = sqrt((p3x - p2x)^2 + (p3y - p2y)^2);

    curvature = 2 * ...
        ((p2x - p1x) * (p3y - p1y) - (p2y - p1y) * (p3x - p1x)) / ...
        (dist1 * dist2 * dist3);

    x_y_yaw_curv_s(idx,3) = yaw;
    x_y_yaw_curv_s(idx,4) = curvature;

    if idx>1
        pcx = x_y_arr(idx, 1);
        pcy = x_y_arr(idx, 2);
        ppx = x_y_arr(idx - 1, 1);
        ppy = x_y_arr(idx - 1, 2);
        
        x_y_yaw_curv_s(idx,5) = x_y_yaw_curv_s(idx-1,5) + ...
            sqrt((pcx - ppx)^2 + (pcy - ppy)^2);
    end
end

x_y_yaw_curv_s(:,1) = x_y_arr(:,1);
x_y_yaw_curv_s(:,2) = x_y_arr(:,2);
end