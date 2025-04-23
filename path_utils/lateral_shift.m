function [x_y_deviated] = lateral_shift(x_y_yaw_arr,deviation)
%lateral_shift Calculate x-y position of the laterally deviated point
%   Given a point that is shifted from the path center-line normally,
%   this function calculates its x-y position
arr_length = length(x_y_yaw_arr(:,1));
x_y_deviated = zeros(arr_length());
for idx = 1:arr_length
    cx = x_y_yaw_arr(idx,1);
    cy = x_y_yaw_arr(idx,2);
    cyaw = x_y_yaw_arr(idx,3);

    dyaw = cyaw + pi/2;
    d = deviation(idx);
    dx = cx + d * cos(dyaw);
    dy = cy + d * sin(dyaw);

    x_y_deviated(idx, 1) = dx;
    x_y_deviated(idx, 2) = dy;
end
end