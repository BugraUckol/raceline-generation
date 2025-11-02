%% Prep
% clc, clear, close all
setDefaultFigureProperties()

%% Curve Generation
% Parameter
num_of_samples = 1000;
t_arr = linspace(0, 2*pi, num_of_samples);  % parametric variable

% Parameters to control shape
a = 3;      % major amplitude (horizontal)
b = 3;    % vertical amplitude
c = 5;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = a * (sin(t_arr) + 2 * sin(2 * t_arr));
y = b * (cos(t_arr) - 2 * cos(2 * t_arr));
z = c * ( - sin(3 * t_arr));
r = [x; y; z];

% Symbolic equations of the curve
t_s = sym("t_s","real");
xs = a * (sin(t_s) + 2 * sin(2 * t_s));
ys = b * (cos(t_s) - 2 * cos(2 * t_s));
zs = c * ( - sin(3 * t_s));
rs = [xs, ys, zs]';

% Symbolic differentiation
r_d = simplify(diff(rs, t_s));

% Creating the tangents array and first normal vector
T_arr = double(subs(r_d/norm(r_d),t_s,t_arr));
nulls = null(T_arr(:,1)');
v1 = nulls(:,1);
V1_arr = [v1];
V2_arr = [cross(T_arr(:,1), v1)];

% Extraction of the arclength
s_arr = [0];
w_arr = [];
for i = 1:num_of_samples-1
    ds = norm(r(:,i+1) - r(:,i));
    s_arr = [s_arr, s_arr(end) + ds];
end

% Creation of the frame
alpha = 2.6817;
figure(2), hold on, daspect([1,1,1])
for i = 1:num_of_samples-1
    ds = s_arr(i+1) - s_arr(i);
    B = cross(T_arr(:,i), T_arr(:,i+1));
    if norm(B) == 0
        V_next = V1_arr(:,end);
        w_dagger_i = [0,0,0]';
        alpha_increment = (alpha / s_arr(end)) * ds;
        spin_rodrigues_vector = T_arr(:,i) * tan(0.5 * alpha_increment);
        rodrigues_vector = spin_rodrigues_vector;
        R_spin = rod2dcm(rodrigues_vector');
        V_next = R_spin' * V_next;
        w_dagger_i = w_dagger_i + alpha_increment * T_arr(:,i+1) / ds;
    else
        B = B / norm(B);
        theta = acos(dot(T_arr(:,i), T_arr(:,i+1)));
        rodrigues_vector = B * tan(0.5 * theta); 
        alpha_increment = (alpha / s_arr(end)) * ds;
        spin_rodrigues_vector = T_arr(:,i) * tan(0.5 * alpha_increment);
        rodrigues_vector = rodrigues_vector + spin_rodrigues_vector;
        R = rod2dcm(rodrigues_vector');
        V_next = R' * V1_arr(:,end);
        w_dagger_i = theta * B / ds + alpha_increment * T_arr(:,i+1) / ds;
    end

    w_dagger_b = [T_arr(:,i+1), V_next, cross(T_arr(:,i+1), V_next)]' * w_dagger_i;

    w_arr = [w_arr, w_dagger_b];
    V1_arr = [V1_arr, V_next];
    V2_arr = [V2_arr, cross(T_arr(:,i+1), V_next)];
    p = double(subs(rs, t_s, t_arr(i+1)));

    if mod(i,5) == 0
        pt = [p, p + T_arr(:,i+1)];
        pn = [p, p + V_next];
        pb = [p, p + cross(T_arr(:,i+1), V_next)];
        plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 2, 'Color', 'r');
        plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 2, 'Color', 'g', 'DisplayName', num2str(i));
        plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 2, 'Color', 'b');
    end
end

w_arr = [w_arr(:,end), w_arr];
figure(1), hold on
plot(s_arr, w_arr)













