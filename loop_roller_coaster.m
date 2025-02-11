clear;
clc;
close all;

%% Loop parameters
h = 75; % height of bottom of loop
radius_loop = 21; % radius of Loop
%% creating graph of loop and transition
x_val = -1;
end_of_tan_line = 50;%adjust this variable later when know more about track
%set to same y value that banked turns ends at
arc_length_loop = 2 * pi * radius_loop;
%when calculating s for loop use this variable instead of integration
height_of_center = h + radius_loop;
figure()
%these plots will be made parametrically
hold on
t1 = 0:pi/16:(2*pi)+(5*pi/4);
x_pos1 = x_val*ones(size(t1));
%parametric values for loop making one full roation and an 1/8 of a roation
lower_bound_t2 = -1*radius_loop*sin(5*pi/4);%finds where to start the tan line
t2 = lower_bound_t2:pi/16:end_of_tan_line;
x_pos2 = x_val*ones(size(t2));
y1 = -1*radius_loop * sin(t1);
z1 = radius_loop * cos(t1);
% y1 and z1 are equations for loop, x should stay constant
y2 = t2;
z2 = t2-radius_loop*sqrt(2);
%y2 and z2 are equations for transitions x should stay constant
set(gca, 'YDir', 'reverse'); % Flips the y-axis
plot3(x_pos1,y1,z1);
plot3(x_pos2,y2,z2);
view(3)
xlabel("x axis");
ylabel("y axis");
zlabel("z axis");
grid on
%% G forces as a function of S plot
% Define vertical G-force function
vertical_g = @(s) (250 - 2 * h) / radius_loop + 3 * cos(s / radius_loop) - 2;
max_vert_g = ((250-2*h)/radius_loop)+1;
min_down_g = ((250-2*h)/radius_loop)-5;
% Create a figure and set up subplots
figure();
hold on

% First subplot: Original plot with vertical G-force
subplot(3,1,1);
fplot(vertical_g, 'LineWidth', 2);
xlim([-10, 2 * pi * radius_loop + 10]);
ylim([-2, 7]);
yline(6, 'r--', 'LineWidth', 2);
yline(-1, 'm--', 'LineWidth', 2);
legend("Vertical G force as a function of S", "Max vertical G", "Max downward G", "Location", "northeast");
text(25, 5.5, "Max Vertical G = 5.76");
text(10, -0.5, "Max Downward G = 0.2381");
title('Vertical G-Force vs Arc Length');
xlabel('Arc Length (m)');
ylabel('G-Force');
%subplot lateral g force
subplot(3,1,2);
yline(0, 'b-', 'LineWidth', 2);
title('Reference Line at 0');
xlabel('Arc Length (m)');
ylabel('G-Force');
yline(-3, 'r--', 'LineWidth', 2);
yline(3, 'm--', 'LineWidth', 2);
xlim([-10, 2 * pi * radius_loop + 10]);
ylim([-4, 4]);
legend("Lateral G force as a function of S", "Max left G", "Max right G", "Location", "northeast");
title('Lateral G-Force vs Arc Length');
%sub plot forward/backward g force
subplot(3,1,3);
yline(0, 'b-', 'LineWidth', 2);
xlabel('Arc Length (m)');
ylabel('G-Force');
yline(5, 'r--', 'LineWidth', 2);
yline(-4, 'm--', 'LineWidth', 2);
xlim([-10, 2 * pi * radius_loop + 10]);
ylim([-6, 6]);
legend("Foward/Backward G force as a function of S", "Max Forward G", "Max Backward G", "Location", "northeast");
title('Forward/Backward G-Force vs Arc Length');

%Zero g parabola
figure();
hold on
subplot(3,1,1);
yline(0,'b-','LineWidth',2);
xlabel("Arc Length (m)");
ylabel("G-Force");
yline(6, 'r--', 'LineWidth', 2);
yline(-1, 'm--', 'LineWidth', 2);
legend("Vertical G force as a function of S", "Max vertical G", "Max downward G", "Location", "northeast");
title('Vertical G-Force vs Arc Length');
%subplot lateral g force
subplot(3,1,2);
yline(0, 'b-', 'LineWidth', 2);
title('Reference Line at 0');
xlabel('Arc Length (m)');
ylabel('G-Force');
yline(-3, 'r--', 'LineWidth', 2);
yline(3, 'm--', 'LineWidth', 2);
ylim([-4, 4]);
legend("Lateral G force as a function of S", "Max left G", "Max right G", "Location", "northeast");
title('Lateral G-Force vs Arc Length');
subplot(3,1,3);
%subplot tengential g force
yline(0, 'b-', 'LineWidth', 2);
xlabel('Arc Length (m)');
ylabel('G-Force');
yline(5, 'r--', 'LineWidth', 2);
yline(-4, 'm--', 'LineWidth', 2);
ylim([-6, 6]);
legend("Foward/Backward G force as a function of S", "Max Forward G", "Max Backward G", "Location", "northeast");
title('Forward/Backward G-Force vs Arc Length');
