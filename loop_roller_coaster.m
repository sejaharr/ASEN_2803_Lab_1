clear;
clc;
close all;

% Define parameters
h = 75; % height of bottom of loop
radius_loop = 21; % radius of Loop
arc_length_loop = 2 * pi * radius_loop;
height_of_center = 75 + radius_loop;
x = @(t) radius_loop * sin(t);
y = @(t) radius_loop * cos(t);

% Define vertical G-force function
vertical_g = @(s) (250 - 2 * h) / radius_loop + 3 * cos(s / radius_loop) - 2;

% Create a figure and set up subplots
figure;
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
xlabel('Arc Length');
ylabel('G-Force');

% Second subplot: Horizontal line at 0
subplot(3,1,2);
yline(0, 'b-', 'LineWidth', 2);
title('Reference Line at 0');
xlabel('Arc Length');
ylabel('G-Force');
xlim([-10, 2 * pi * radius_loop + 10]);
ylim([-1, 1]);
title('Lateral G-Force vs Arc Length');

% Third subplot: Another horizontal line at 0
subplot(3,1,3);
yline(0, 'b-', 'LineWidth', 2);
title('Another Reference Line at 0');
xlabel('Arc Length');
ylabel('G-Force');
xlim([-10, 2 * pi * radius_loop + 10]);
ylim([-1, 1]);
title('Forward/Backward G-Force vs Arc Length');


hold off;


