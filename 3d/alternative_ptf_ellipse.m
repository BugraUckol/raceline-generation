%% Prep
% clc, clear, close all
setDefaultFigureProperties()

%% Curve Generation
% Parameter
num_of_samples = 500;
t_arr = linspace(0, 2*pi, num_of_samples);  % parametric variable

% Parameters to control shape
a = 5;      % major amplitude (horizontal)
b = 2;    % vertical amplitude
c = 0;    % depth amplitude (3D deviation)

% Parametric equations for 3D figure-eight curve
x = a * cos(t_arr);
y = b * sin(t_arr);
z = 0 * t_arr;
r = [x; y; z];

% Symbolic equations of the curve
t_s = sym("t_s","real");
xs = a * cos(t_s);
ys = b * sin(t_s);
zs = c * t_s;
rs = [xs, ys, zs]';

% Symbolic differentiation
r_d = simplify(diff(rs, t_s));

% Creating the tangents array and first normal vector
T_arr = double(subs(r_d/norm(r_d),t_s,t_arr));
nulls = null(T_arr(:,1)');
v1 = nulls(:,1);
V1_arr = [v1];
V2_arr = [cross(T_arr(:,1), v1)];
[yaw, pitch, roll] = ...
    dcm2angle( [T_arr(:,1),v1,V2_arr(:,1)]', 'zyx', 'robust');
yaw_arr = [yaw];
pitch_arr = [pitch];
roll_arr = [roll];


% Extraction of the arclength
s_arr = [0];
w_arr = [];
for i = 1:num_of_samples-1
    ds = norm(r(:,i+1) - r(:,i));
    s_arr = [s_arr, s_arr(end) + ds];
end

% Creating circle edge around the tangent vector
radius = 0.65;
circle_edges = radius * [0*cos(-0.1:0.1:2*pi); cos(-0.1:0.1:2*pi); sin(-0.1:0.1:2*pi)];
circle_edge_arr = r(:,1) + [T_arr(:,1), v1, V2_arr(:,1)] * circle_edges;

% Creation of the frame
alpha = 0;
figure(2), hold on, daspect([1,1,1])
h(1) = plot3(x, y, z, 'k', 'LineWidth', 2, 'DisplayName','$$\mathcal{O}_{Reference}$$');
axis equal;
xlabel('East'); ylabel('North'); zlabel('Up');
title('3D Tube');
view(75, 15); % adjust view angle for clarity
hold on
box on;
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

    if mod(i,15) == 0
        pt = [p, p + T_arr(:,i+1)*0.5];
        pn = [p, p + V_next*0.5];
        pb = [p, p + cross(T_arr(:,i+1), V_next)*0.5];
        plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 2, 'Color', 'r');
        plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 2, 'Color', 'g', 'DisplayName', num2str(i));
        plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 2, 'Color', 'b');

        circle_edges_i = p + [T_arr(:,i+1), V_next, cross(T_arr(:,i+1), V_next)] * circle_edges;
        circle_edge_arr = cat(3,circle_edge_arr,circle_edges_i);
        plot3(circle_edges_i(1,:),circle_edges_i(2,:),circle_edges_i(3,:),'k','LineWidth',0.2);
    end

    if i == 1
        pt = [p, p + T_arr(:,i+1)*0.5];
        pn = [p, p + V_next*0.5];
        pb = [p, p + cross(T_arr(:,i+1), V_next)*0.5];
        circle_edges_i = p + [T_arr(:,i+1), V_next, cross(T_arr(:,i+1), V_next)] * circle_edges;
        h(2) = plot3(circle_edges_i(1,:),circle_edges_i(2,:),circle_edges_i(3,:),'c','LineWidth',4,'DisplayName','Start');
        h(4) = plot3(pt(1,:), pt(2,:), pt(3,:), 'LineWidth', 2, 'Color', 'r','DisplayName','$$\hat t$$');
        h(5) = plot3(pn(1,:), pn(2,:), pn(3,:), 'LineWidth', 2, 'Color', 'g', 'DisplayName', num2str(i),'DisplayName','$$\hat{n}_1$$');
        h(6) = plot3(pb(1,:), pb(2,:), pb(3,:), 'LineWidth', 2, 'Color', 'b','DisplayName','$$\hat{n}_2$$');
    elseif i == num_of_samples-1
        circle_edges_i = p + [T_arr(:,i+1), V_next, cross(T_arr(:,i+1), V_next)] * circle_edges;
        h(3) = plot3(circle_edges_i(1,:),circle_edges_i(2,:),circle_edges_i(3,:),'m','LineWidth',4,'DisplayName','Finish')
    end

    [yaw, pitch, roll] = ...
        dcm2angle( [T_arr(:,i+1),V_next,cross(T_arr(:,i+1), V_next)]', 'zyx', 'robust');
    yaw_arr = [yaw_arr, yaw];
    pitch_arr = [pitch_arr, pitch];
    roll_arr = [roll_arr, roll];
end
circle_edge_arr = cat(3,circle_edge_arr,circle_edge_arr(:,:,1));
s = surf(squeeze(circle_edge_arr(1,:,:)), ...
    squeeze(circle_edge_arr(2,:,:)), ...
    squeeze(circle_edge_arr(3,:,:)));
s.FaceAlpha = 0.4;
s.EdgeColor = 'none';
s.FaceColor = [0.5, 0.5, 0.5];
% legend({'Centerline','Start','$$\hat t$$','$$\hat n_1$$','$$\hat n_2$$'},'interpreter','latex')
legend(h,'interpreter','latex')

kappa = [w_arr(:,end), w_arr];
figure(1), hold on
plot(s_arr, kappa)
xlabel('Curvilinear Distance [m]')
ylabel('Value [1/m]')
legend('$$\tau$$','$$\kappa_1$$','$$\kappa_2$$','interpreter','latex')
title('Arclength -vs- Generalized Curvature Profiles of The Curve')
grid minor

x_arr = x;
y_arr = y;
z_arr = z;
tau_arr = w_arr(1,:);
kap1_arr = w_arr(2,:);
kap2_arr = w_arr(3,:);
tt_arr = T_arr;
nn_arr = V1_arr;
bb_arr = V2_arr;

save alternative_ellipse_ptf s_arr x_arr y_arr z_arr tau_arr kap1_arr kap2_arr tt_arr ...
    nn_arr bb_arr yaw_arr pitch_arr roll_arr circle_edge_arr













