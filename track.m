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
theta_f_parab = pi/4;
t_final_parab =3.445;%got this value experimentally
g = 9.81;
arc_length_segments = zeros(7);
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
plot3(x,y,z);

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
plot3(x_trans_1,y_trans_1,z_trans_1);

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
plot3(x_loop,y_loop,z_loop);

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
plot3(x_trans_2,y_trans_2,z_trans_2);

%% Arc Length Transition 2

dx2 = diff(x_trans_2);
dz2 = diff(z_trans_2);
dy2 = diff(y_trans_2);
ds2 = sqrt(dx2.^2 + dy2.^2 + dz2.^2); % Small segment lengths
arc_length_segments(4) = sum(ds2); % Sum of segment lengths

%% Graphing zero g parabola

t_parab = 0:0.1:t_final_parab;
end_trans_2.x = end_loop.x;
end_trans_2.z = trans_2_end-radius_loop*sqrt(2)+h_center;
end_trans_2.y = -1*trans_2_end+end_trans_1.y;
first_term = tan(theta_i_parab)*t_parab;
v_o_parab = sqrt(2*g*(initial_height-end_trans_2.z));
v_o_z = v_o_parab*sin(theta_i_parab);
v_o_y = v_o_parab*cos(theta_f_parab);
z_parab = v_o_z*t_parab-0.5*g*t_parab.^2+end_trans_2.z;
x_parab = end_trans_2.x*ones(size(t_parab));
y_parab = -1*v_o_y*t_parab+end_trans_2.y;
plot3(x_parab,y_parab,z_parab);

%% Arc Length Zero G Parabola 

dx_parab = diff(x_parab);
dz_parab = diff(z_parab);
dy_parab = diff(y_parab);
ds_parab = sqrt(dx_parab.^2 + dy_parab.^2 + dz_parab.^2); % Small segment lengths
arc_length_segments(5) = sum(ds_parab); % Sum of segment lengths

%% Graphing Transition 3
%{
end_parab.z = v_o_z*(t_final_parab)-0.5*g*(t_final_parab).^2+end_trans_2.z;
end_parab.x = max(x_parab);
end_parab.y = -1*v_o_y*t_final_parab+end_trans_2.y;
radius_trans_3 = abs(end_parab.z)*4;
t_trans_3 = 7*pi/4:pi/64:2*pi;
y_trans_3 = -1*(radius_trans_3*sin(t_trans_3))+end_parab.y;
z_trans_3 = -1*(radius_trans_3*cos(t_trans_3))+radius_trans_3;
x_trans_3 = end_parab.x*ones(size(t_trans_3));
plot3(x_trans_3,y_trans_3,z_trans_3);
%}
%% Finding Final Arc Length
arc_length = max(sum(arc_length_segments));

%% Graphing Essentials
xlabel("x axis");
ylabel("y axis");
zlabel("z axis");
view(3)
grid on
