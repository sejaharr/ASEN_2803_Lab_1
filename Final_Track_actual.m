clear; clc; close all;
figure();
hold on;

%% Constants
radius_helix = 30;         % 20; 15;
n = 2.5;                   % 3.5; 4.5;  % number of loops
initial_height = 125;
c = 1;                     % 1.5;  % rate at which helix descends
radius_loop = 7;           % 10; 16; 20;
trans_1_end = 40;          % 35;
trans_2_end = 3*pi;
lower_bound_t_trans_2 = -radius_loop*sin(5*pi/4);
theta_i_parab = pi/4;
t_final_parab = 4.9;       % 5.5; % got this value experimentally
g = 9.81;
arc_length_segments = zeros(8,1);
banked_turn_2 = 0;         % 1 means do the optional banked turn
radius_opt_banked_turn = 50;
length_of_breaking_section = 150;

%% Debug Plane at z = 125
%{
[x_debug, y_debug] = meshgrid(-100:5:100, -100:5:100);
z_debug = 125*ones(size(x_debug));
surf(x_debug, y_debug, z_debug, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%}

%% Graphing Helix
t_helix = 0:pi/32:(2*pi)*n;
x = radius_helix*cos(t_helix);
y = radius_helix*sin(t_helix);
z = -c*t_helix + initial_height;
plotSegment3D(x, y, z, initial_height, g);

% Calculating lateral G for helix
V = sqrt(2 * g * (initial_height - z));
lateral_g_helix = V.^2 / (radius_helix * g);
% vertical G is constant = 1

% Computing Arc Length of Helix
arc_length_segments(1) = computeArcLength(x, y, z);

%% Graphing Transition 1
end_helix.x = radius_helix*cos(2*pi*n);
end_helix.y = radius_helix*sin(2*pi*n);
end_helix.z = -c*(2*pi*n) + initial_height;
t_trans_1 = 0:1:trans_1_end;
x_trans_1 = end_helix.x * ones(size(t_trans_1));
z_trans_1 = end_helix.z * ones(size(t_trans_1));
y_trans_1 = -t_trans_1 + end_helix.y;
plotSegment3D(x_trans_1, y_trans_1, z_trans_1, initial_height, g);
arc_length_segments(2) = abs(trans_1_end);

%% Graphing Loop
end_trans_1.x = max(x_trans_1);
end_trans_1.z = max(z_trans_1);
end_trans_1.y = -trans_1_end + end_helix.y;
h_b = end_trans_1.z;
h_center = h_b + radius_loop;
t_loop = 0:pi/32:(2*pi)+(5*pi/4);
x_loop = end_trans_1.x * ones(size(t_loop));
y_loop = -radius_loop*sin(t_loop) + end_trans_1.y;
z_loop = -radius_loop*cos(t_loop) + h_center;
plotSegment3D(x_loop, y_loop, z_loop, initial_height, g);
arc_length_segments(3) = 2*pi*radius_loop;

%% Graphing Transition 2
end_loop.x = max(x_loop);
end_loop.z = -radius_loop*cos((2*pi)+(5*pi/4)) + h_center;
end_loop.y = -radius_loop*sin((2*pi)+(5*pi/4)) + end_trans_1.y;
t_trans_2 = lower_bound_t_trans_2:pi/32:trans_2_end;
x_trans_2 = end_loop.x * ones(size(t_trans_2));
y_trans_2 = -t_trans_2 + end_trans_1.y;
z_trans_2 = t_trans_2 - radius_loop*sqrt(2) + h_center;
plotSegment3D(x_trans_2, y_trans_2, z_trans_2, initial_height, g);
arc_length_segments(4) = computeArcLength(x_trans_2, y_trans_2, z_trans_2);

%% Graphing Zero-G Parabola
t_parab = 0:0.001:t_final_parab;
end_trans_2.x = end_loop.x;
end_trans_2.z = trans_2_end - radius_loop*sqrt(2) + h_center;
end_trans_2.y = -trans_2_end + end_trans_1.y;
v_o_parab = sqrt(2*g*(initial_height - end_trans_2.z));
v_o_z = v_o_parab*sin(theta_i_parab);
v_o_y = v_o_parab*cos(theta_i_parab);
z_parab = v_o_z*t_parab - 0.5*g*t_parab.^2 + end_trans_2.z;
x_parab = end_trans_2.x * ones(size(t_parab));
y_parab = -v_o_y*t_parab + end_trans_2.y;
plotSegment3D(x_parab, y_parab, z_parab, initial_height, g);
arc_length_segments(5) = computeArcLength(x_parab, y_parab, z_parab);

%% Graphing Transition 3
end_parab.z = v_o_z*t_final_parab - 0.5*g*t_final_parab^2 + end_trans_2.z;
end_parab.x = max(x_parab);
end_parab.y = -v_o_y*t_final_parab + end_trans_2.y;
d_f_parab = (v_o_y - g*t_final_parab)/v_o_z;
theta_f_parab = atan(abs(d_f_parab));
radius_trans_3 = -end_parab.z/(cos(theta_f_parab)-1);
delta_y = radius_trans_3*sin(theta_f_parab);
t_trans_3 = 0:0.005:theta_f_parab;
x_trans_3 = end_parab.x * ones(size(t_trans_3));
y_trans_3 = radius_trans_3*sin(t_trans_3) + end_parab.y - delta_y;
z_trans_3 = -radius_trans_3*cos(t_trans_3) + (radius_trans_3+1);
plotSegment3D(x_trans_3, y_trans_3, z_trans_3, initial_height, g);
arc_length_segments(6) = radius_trans_3 * theta_f_parab;

%% Graphing Optional Banked Turn
end_trans_3.z = -radius_trans_3*cos(0) + radius_trans_3;
end_trans_3.y = radius_trans_3*sin(0) + end_parab.y - delta_y;
end_trans_3.x = max(x_trans_3);
if banked_turn_2 == 1
    t_opt_banked_turn = pi/2:pi/256:3*pi/2;
    z_opt_banked_turn = end_trans_3.z * ones(size(t_opt_banked_turn));
    x_opt_banked_turn = radius_opt_banked_turn*sin(t_opt_banked_turn) + end_trans_3.x - radius_opt_banked_turn;
    y_opt_banked_turn = radius_opt_banked_turn*cos(t_opt_banked_turn) + end_trans_3.y;
    plotSegment3D(x_opt_banked_turn, y_opt_banked_turn, z_opt_banked_turn, initial_height, g);
    arc_length_segments(7) = pi * radius_opt_banked_turn;
end

%% Graphing Braking Section
velocity_at_start = sqrt(2*g*initial_height);
acceleration = -velocity_at_start/length_of_breaking_section;
t_braking_section = 0:length_of_breaking_section;
if banked_turn_2 == 1
    z_pos_braking = end_trans_3.z * ones(size(t_braking_section));
    x_pos_braking = (radius_opt_banked_turn*sin(3*pi/2) + end_trans_3.x - radius_opt_banked_turn)*ones(size(t_braking_section));
    y_pos_braking = radius_opt_banked_turn*cos(3*pi/2) + end_trans_3.y + t_braking_section;
else
    z_pos_braking = end_trans_3.z * ones(size(t_braking_section));
    x_pos_braking = end_trans_3.x * ones(size(t_braking_section));
    y_pos_braking = end_trans_3.y - t_braking_section;
end

V = velocity_at_start + acceleration*t_braking_section;
scatter3(x_pos_braking, y_pos_braking, z_pos_braking, 20, V, 'filled'); 
arc_length_segments(8) = length_of_breaking_section;

%% Final Arc Length
arc_length = sum(arc_length_segments);

%% Graphing Set Up
xlabel("x axis");
ylabel("y axis");
zlabel("z axis");
view(3);
grid on;
axis equal;

%% G(s) Plotting
% Create a structure array to hold each segment's G-force data.
gData = struct();

% Helix
track_lengths_helix = linspace(0, arc_length_segments(1), length(lateral_g_helix));
gData(1).arc = track_lengths_helix;
gData(1).vertical = ones(size(track_lengths_helix));
gData(1).lateral = lateral_g_helix;
gData(1).tangential = zeros(size(track_lengths_helix));
gData(1).title = 'Helix';

% Transition 1
track_lengths_trans_1 = track_lengths_helix(end):0.1:(track_lengths_helix(end)+arc_length_segments(2));
gData(2).arc = track_lengths_trans_1;
gData(2).vertical = ones(size(track_lengths_trans_1));
gData(2).lateral = zeros(size(track_lengths_trans_1));
gData(2).tangential = zeros(size(track_lengths_trans_1));
gData(2).title = 'Transition 1';

% Loop (using computed values)
track_lengths_loop = track_lengths_trans_1(end):0.1:(track_lengths_trans_1(end)+arc_length_segments(3));
max_vertical = ((250-2*end_trans_1.z)/radius_loop) + 1;
vertical_g_loop = ((250-2*end_trans_1.z)/radius_loop) + 3*cos((track_lengths_loop - track_lengths_trans_1(end))/radius_loop) - 2;
lateral_g_loop = zeros(size(track_lengths_loop));
fb_g_loop = zeros(size(track_lengths_loop));
gData(3).arc = track_lengths_loop;
gData(3).vertical = vertical_g_loop;
gData(3).lateral = lateral_g_loop;
gData(3).tangential = fb_g_loop;
gData(3).title = 'Loop';

% Transition 2
track_lengths_trans_2 = track_lengths_loop(end):0.1:(track_lengths_loop(end)+arc_length_segments(4));
gData(4).arc = track_lengths_trans_2;
gData(4).vertical = ones(size(track_lengths_trans_2));
gData(4).lateral = zeros(size(track_lengths_trans_2));
gData(4).tangential = -ones(size(track_lengths_trans_2));
gData(4).title = 'Transition 2';

% Parabola
track_lengths_parab = track_lengths_trans_2(end):0.1:(track_lengths_trans_2(end)+arc_length_segments(5));
gData(5).arc = track_lengths_parab;
gData(5).vertical = zeros(size(track_lengths_parab));
gData(5).lateral = zeros(size(track_lengths_parab));
gData(5).tangential = zeros(size(track_lengths_parab));
gData(5).title = 'Parabola';

% Transition 3
track_lengths_trans_3 = track_lengths_parab(end):0.1:(track_lengths_parab(end)+arc_length_segments(6));
theta_f_parab = 2*pi-theta_f_parab;
vertical_g_trans_3 = ((250-end_trans_3.z)/radius_trans_3) + 3*cos(((theta_f_parab*radius_trans_3)+(track_lengths_trans_3 - track_lengths_parab(end)))/radius_trans_3) - 2;
lateral_g_trans_3 = zeros(size(track_lengths_trans_3));
fb_g_trans_3 = zeros(size(track_lengths_trans_3));
gData(6).arc = track_lengths_trans_3;
gData(6).vertical = vertical_g_trans_3;
gData(6).lateral = lateral_g_trans_3;
gData(6).tangential = fb_g_trans_3;
gData(6).title = 'Transition 3';

% Braking Section
track_lengths_brake = track_lengths_trans_3(end):0.1:(track_lengths_trans_3(end)+arc_length_segments(8));
vertical_g_brake = ones(size(track_lengths_brake));
lateral_g_brake = zeros(size(track_lengths_brake));
accel_brake = (velocity_at_start^2) / (2*arc_length_segments(8));
coeff_friction = 0.4;
fb_g_brake = (accel_brake / (coeff_friction*g)) * ones(size(track_lengths_brake));
gData(7).arc = track_lengths_brake;
gData(7).vertical = vertical_g_brake;
gData(7).lateral = lateral_g_brake;
gData(7).tangential = fb_g_brake;
gData(7).title = 'Braking Section';

% Entire Track data (concatenated)
total_track_length = [gData(1).arc, gData(2).arc, gData(3).arc, gData(4).arc, gData(5).arc, gData(6).arc, gData(7).arc];
total_vertical_g = [gData(1).vertical, gData(2).vertical, gData(3).vertical, gData(4).vertical, gData(5).vertical, gData(6).vertical, gData(7).vertical];
total_lateral_g = [gData(1).lateral, gData(2).lateral, gData(3).lateral, gData(4).lateral, gData(5).lateral, gData(6).lateral, gData(7).lateral];
total_fb_g = [gData(1).tangential, gData(2).tangential, gData(3).tangential, gData(4).tangential, gData(5).tangential, gData(6).tangential, gData(7).tangential];
cumulative_arc_lengths = cumsum(arc_length_segments);

% Plot Entire Track G-forces
plotGForceSegment(total_track_length, total_vertical_g, total_lateral_g, total_fb_g, cumulative_arc_lengths, "Entire Track");

% Plot individual segment G-forces in separate figures
for k = 1:length(gData)
    plotGForceSegment(gData(k).arc, gData(k).vertical, gData(k).lateral, gData(k).tangential, [], gData(k).title);
end

%% Local Functions

function plotSegment3D(x, y, z, initHeight, g)
    % Compute velocity (based on difference from initHeight) for coloring
    V = sqrt(2 * g * (initHeight - z));
    scatter3(x, y, z, 20, V, 'filled');
    colormap(jet);
    colorbar;
end

function L = computeArcLength(x, y, z)
    dx = diff(x);
    dy = diff(y);
    dz = diff(z);
    ds = sqrt(dx.^2 + dy.^2 + dz.^2);
    L = sum(ds);
end

function plotGForceSegment(arc, vertical, lateral, tangential, segBoundaries, segTitle)
    figure();
    % Vertical G-force subplot
    subplot(3,1,1);
    plot(arc, vertical, 'b-', 'LineWidth', 2);
    ylim([-2, 7]);
    yline(6, 'r--', 'LineWidth', 2);
    yline(-1, 'm--', 'LineWidth', 2);
    xlabel('Arc Length (m)');
    ylabel('G-Force');
    title(['Vertical G-Force: ' segTitle]);
    grid on;
    
    % Lateral G-force subplot
    subplot(3,1,2);
    plot(arc, lateral, 'b-', 'LineWidth', 2);
    ylim([-4, 4]);
    yline(-3, 'r--', 'LineWidth', 2);
    yline(3, 'm--', 'LineWidth', 2);
    xlabel('Arc Length (m)');
    ylabel('G-Force');
    title(['Lateral G-Force: ' segTitle]);
    grid on;
    
    % Tangential G-force subplot
    subplot(3,1,3);
    plot(arc, tangential, 'b-', 'LineWidth', 2);
    ylim([-5, 6]);
    yline(5, 'r--', 'LineWidth', 2);
    yline(-4, 'm--', 'LineWidth', 2);
    xlabel('Arc Length (m)');
    ylabel('G-Force');
    title(['Tangential G-Force: ' segTitle]);
    grid on;
    
    % draws vertical lines on all subplots
    if ~isempty(segBoundaries)
        for i = 1:length(segBoundaries)
            subplot(3,1,1);
            xline(segBoundaries(i), 'k--', 'LineWidth', 1);
            subplot(3,1,2);
            xline(segBoundaries(i), 'k--', 'LineWidth', 1);
            subplot(3,1,3);
            xline(segBoundaries(i), 'k--', 'LineWidth', 1);
        end
    end
    sgtitle(['G-forces for ', segTitle]);
end
