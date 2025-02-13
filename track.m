clear;clc;close all;
figure()
hold on
%% Constants

radius_helix = 10;
n = 5.5;
initial_height = 125;
c=1.5;
radius_loop = 21;
trans_1_end = 30;
trans_2_end = 10*pi;
lower_bound_t_trans_2 = -1*radius_loop*sin(5*pi/4);
theta_i_parab = pi/4;
t_final_parab =5.5;%got this value experimentally
g = 9.81;
arc_length_segments = zeros(8);
banked_turn_2 = 1;% one means do the optional banked turn
radius_opt_banked_turn = 50;
length_of_breaking_section = 150;


%% PLane at z =125 for debugging 
%{
% Define x and y ranges
[x_debug, y_debug] = meshgrid(-100:5:100, -100:5:100); % Adjust the range as needed

% Set constant z-value for the plane
z_debug = 125 * ones(size(x_debug));

% Plot the plane
surf(x_debug, y_debug, z_debug, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Semi-transparent plane
%}
%% Graphing Helix

t_helix = 0:pi/32:(2*pi)*n;
x=radius_helix*cos(t_helix);
y=radius_helix*sin(t_helix);
z=-c*t_helix+initial_height;
V = sqrt(2 * g * (125 - z));
scatter3(x, y, z, 20, V, 'filled'); 
colormap(jet);
colorbar; 
 
%% Computing Arc Length of Helix
dx = diff(x);
dy = diff(y);
dz = diff(z);

ds = sqrt(dx.^2 + dy.^2 + dz.^2); % Small segment lengths
arc_length_segments(1) = sum(ds); % Sum of segment lengths

%% Graphing Transition 1

end_helix.y = radius_helix*sin((2*pi)*n);
end_helix.x = radius_helix*cos((2*pi)*n);
end_helix.z = -c*(2*pi*n)+initial_height;
t_trans_1 = 0:1:trans_1_end;
x_trans_1 = end_helix.x*ones(size(t_trans_1));
z_trans_1 = end_helix.z*ones(size(t_trans_1));
y_trans_1 = -1*t_trans_1+end_helix.y;
V = sqrt(2 * g * (125 - z_trans_1));
scatter3(x_trans_1, y_trans_1, z_trans_1, 20, V, 'filled'); 
colormap(jet);
colorbar; 

%% Arc Length Transition 1

arc_length_segments(2) = abs(trans_1_end);

%% Graphing Loop

end_trans_1.x = max(x_trans_1);
end_trans_1.z = max(z_trans_1);
end_trans_1.y = -1*trans_1_end+end_helix.y;
h_b = end_trans_1.z;
h_center = h_b+radius_loop;
t_loop = 0:pi/32:(2*pi)+(5*pi/4);
x_loop = end_trans_1.x*ones(size(t_loop));
y_loop = -1*radius_loop*sin(t_loop)+end_trans_1.y;
z_loop = -1*radius_loop*cos(t_loop)+h_center;
V = sqrt(2 * g * (125 - z_loop));
scatter3(x_loop, y_loop, z_loop, 20, V, 'filled'); 
colormap(jet);
colorbar; 

%% Arc Length Loop

arc_length_segments(3) = 2*pi*radius_loop;

%% Graphing Transition 2

end_loop.x = max(x_loop);
end_loop.z = -1*radius_loop*cos((2*pi)+(5*pi/4))+h_center;
end_loop.y = -1*radius_loop*sin((2*pi)+(5*pi/4))+end_trans_1.y;
t_trans_2 = lower_bound_t_trans_2:pi/32:trans_2_end;
x_trans_2 = end_loop.x*ones(size(t_trans_2));
y_trans_2 = -1*t_trans_2+end_trans_1.y;
z_trans_2 = t_trans_2-radius_loop*sqrt(2)+h_center;
V = sqrt(2 * g * (125 - z_trans_2));
scatter3(x_trans_2, y_trans_2, z_trans_2, 20, V, 'filled'); 
colormap(jet);
colorbar; 

%% Arc Length Transition 2

dx2 = diff(x_trans_2);
dz2 = diff(z_trans_2);
dy2 = diff(y_trans_2);
ds2 = sqrt(dx2.^2 + dy2.^2 + dz2.^2); 
arc_length_segments(4) = sum(ds2); 

%% Graphing zero g parabola

t_parab = 0:0.001:t_final_parab;
end_trans_2.x = end_loop.x;
end_trans_2.z = trans_2_end-radius_loop*sqrt(2)+h_center;
end_trans_2.y = -1*trans_2_end+end_trans_1.y;
first_term = tan(theta_i_parab)*t_parab;
v_o_parab = sqrt(2*g*(initial_height-end_trans_2.z));
v_o_z = v_o_parab*sin(theta_i_parab);
v_o_y = v_o_parab*cos(theta_i_parab);
z_parab = v_o_z*t_parab-0.5*g*t_parab.^2+end_trans_2.z;
x_parab = end_trans_2.x*ones(size(t_parab));
y_parab = -1*v_o_y*t_parab+end_trans_2.y;
V = sqrt(2 * g * (125 - z_parab));
scatter3(x_parab, y_parab, z_parab, 20, V, 'filled'); 
colormap(jet);
colorbar; 

%% Arc Length Zero G Parabola 

dx_parab = diff(x_parab);
dz_parab = diff(z_parab);
dy_parab = diff(y_parab);
ds_parab = sqrt(dx_parab.^2 + dy_parab.^2 + dz_parab.^2); % Small segment lengths
arc_length_segments(5) = sum(ds_parab); % Sum of segment lengths

%% Graphing Transition 3

end_parab.z = v_o_z*(t_final_parab)-0.5*g*(t_final_parab).^2+end_trans_2.z;
end_parab.x = max(x_parab);
end_parab.y = -1*v_o_y*t_final_parab+end_trans_2.y;
d_f_parab = (v_o_y-g*(t_final_parab))/v_o_z;
theta_f_parab = atan(abs(d_f_parab));
radius_trans_3 = (-1*end_parab.z)/(cos(theta_f_parab)-1);
delta_y = radius_trans_3*sin(theta_f_parab);
t_trans_3 = 0:0.005:theta_f_parab;
x_trans_3 = end_parab.x*ones(size(t_trans_3));
y_trans_3 = radius_trans_3*sin(t_trans_3)+end_parab.y-delta_y;
z_trans_3 = -1*radius_trans_3*cos(t_trans_3)+(radius_trans_3+1);
V = sqrt(2 * g * (125 - z_trans_3));
scatter3(x_trans_3, y_trans_3, z_trans_3, 20, V, 'filled'); 
colormap(jet);
colorbar; 

%% Arc Length Transition 3

arc_length_segments(6) = pi*radius_trans_3*(theta_f_parab);

%% Graphing Optional Banked Turn

if banked_turn_2==1
    end_trans_3.z = -1*radius_trans_3*cos(0)+(radius_trans_3);
    end_trans_3.y = radius_trans_3*sin(0)+end_parab.y-delta_y;
    end_trans_3.x = max(x_trans_3);
    t_opt_banked_turn = pi/2:pi/256:3*pi/2;
    z_opt_banked_turn = end_trans_3.z*ones(size(t_opt_banked_turn));
    x_opt_banked_turn = radius_opt_banked_turn*sin(t_opt_banked_turn)+end_trans_3.x-radius_opt_banked_turn;
    y_opt_banked_turn = radius_opt_banked_turn*cos(t_opt_banked_turn)+end_trans_3.y;
    V = sqrt(2 * g * (125 - z_opt_banked_turn));
    scatter3(x_opt_banked_turn, y_opt_banked_turn, z_opt_banked_turn, 20, V, 'filled'); 
    colormap(jet);
    colorbar;
    arc_length_segments(7) = pi*radius_opt_banked_turn;
end

%% Graphing Braking Section

velocity_at_start = sqrt(2*g*initial_height);
acceleration = -1*velocity_at_start/length_of_breaking_section;
t_braking_section = 0:length_of_breaking_section;
if banked_turn_2 ==1
    z_pos_braking = end_trans_3.z*ones(size(t_braking_section));
    x_pos_braking = (radius_opt_banked_turn*sin(3*pi/2)+end_trans_3.x-radius_opt_banked_turn)*ones(size(t_braking_section));
    y_pos_braking = radius_opt_banked_turn*cos(3*pi/2)+end_trans_3.y+t_braking_section;
end
V = velocity_at_start+acceleration*t_braking_section;
v_end = velocity_at_start+acceleration*100;
scatter3(x_pos_braking, y_pos_braking, z_pos_braking, 20, V, 'filled'); 
colormap(jet);
colorbar;

%% Finding Arc Lenght Braking Section

arc_length_segments(8) = length_of_breaking_section;
%% Finding Final Arc Length
arc_length = max(sum(arc_length_segments));

%% Graphing Essentials
xlabel("x axis");
ylabel("y axis");
zlabel("z axis");
view(3)
grid on
axis equal