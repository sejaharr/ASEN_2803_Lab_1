%% Composite Track Simulation: G Force vs. Distance
clear; 
clc; 
close all;

g = 9.81;          % gravitational acceleration (m/s^2)
h_init = 125;      % initial height (m)
n = 2000;          % number of sample points

% Section time boundaries 
% change these values to space out sections

T1     = 3*pi;     % <-- End of helix (increase T1 to extend helix)
Delta1 = 3;        % <-- Duration of helix→loop transition (increase for longer transition)
T2     = T1 + Delta1;  
T_loop = 2*pi;     % <-- Duration of loop (increase for a longer loop)
T3     = T2 + T_loop;  
Delta2 = 5;        % <-- Duration of loop→projectile transition (increase for longer transition)
T3p    = T3 + Delta2;
T_parab = 4;       % <-- Duration of projectile (parabolic) segment
T4     = T3p + T_parab;
Delta3 = 4;        % Transition duration from parabola to braking
T4p = T4 + Delta3; % Start of the braking section
T_brake = 5;       % Duration of the braking section
T5 = T4p + T_brake; % End of braking section

BrakeOffset = 100;

theta0 = pi/6;     % Starting loop parameter (radians)
k_loop = 1;        % Loop parameter rate



% Create simulation time vector
t_values = linspace(0, T5, n);

% Composite track function (scalar input)
compositeTrack = @(t) compTrack(t, T1, T2, T3, T3p, T4, T4p, T5, theta0, k_loop, g, h_init);

% Evaluate the track at each t (r_all is 3 x n)
r_cell = arrayfun(@(tt) compositeTrack(tt), t_values, 'UniformOutput', false);
r_all = cell2mat(r_cell);
x_t = r_all(1,:);
y_t = r_all(2,:);
z_t = r_all(3,:);

% bank angle (in degrees)
beta_values = 10 * ones(size(t_values));

% Compute track arc length, G force, and speed
[s_track, G_force, ~, ~, ~, v_t, ~, ~] = computeGForces(compositeTrack, h_init, n, t_values, beta_values);

%Calculating G-forces
Gforce = gforce(t_values, T1, T2,T3, T3p, T4, theta0, k_loop, v_t);

% Plot G Force vs. Distance
figure;
plot(s_track, G_force, 'b-', 'LineWidth', 2);
xlabel('Distance (m)'); ylabel('Normal Force (G)');
title('G Force vs. Distance');
grid on;

% Plot 3D track for visual reference
figure;
plot3(x_t, y_t, z_t, 'k-', 'LineWidth', 2);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Composite 3D Track');
axis equal; grid on; rotate3d on;

%Plotting Gforce
figure;
plot(s_track,Gforce)

% compTrack: composite track function with two transitions
function r = compTrack(t, T1, T2, T3, T3p, T4, T4p, ~, theta0, k_loop, g, h_init)
    theta_end = theta0 + (T3 - T2) * k_loop;
    r4 = [10 * sin(theta_end) + 35; 0; h_init - 10 * cos(theta_end)] + [20; 0; 15];
    v4 = [10; 0; 15];
    t_par_final = T4 - T3p;
    r4_end = r4 + v4 * t_par_final + [0; 0; -0.5 * g * (t_par_final^2)];
    v4_end = v4 + [0; 0; -g * t_par_final];
    
    if t < T1
        % Segment 1: Helix (banked turn)
        r = [12*cos(t); 12*sin(t); h_init - (5/(2*pi))*t];
    
    elseif t < T2
        
        % Segment 2: Transition from helix to loop
        L = T2 - T1;
        s_val = (t - T1) / L;
        
        % Helix endpoint at T1
        r1 = [12*cos(T1); 12*sin(T1); h_init - (5/(2*pi))*T1];
        r1p = [-12*sin(T1); 12*cos(T1); -5/(2*pi)];
        
        % Loop start target at theta0
        r2 = [10*sin(theta0)+25; 0; h_init - 10*cos(theta0)];
        r2p_theta = [10*cos(theta0); 0; 10*sin(theta0)];
        r2p = k_loop * r2p_theta;
        a0 = r1; 
        a1 = L * r1p;
        D = r2 - r1 - a1;
        a3 = L * r2p - a1 - 2*D;
        a2 = D - a3;
        r = a0 + a1*s_val + a2*s_val^2 + a3*s_val^3;
    
    elseif t < T3
        % Segment 3: Loop
        theta = theta0 + (t - T2)*k_loop;
        r = [10*sin(theta)+25; 0; h_init - 10*cos(theta)];
    
    elseif t < T3p
        % Segment 4: Transition from loop to projectile
        L2 = T3p - T3;
        s_val2 = (t - T3) / L2;
        
        % Loop endpoint at T3
        theta_end = theta0 + (T3 - T2)*k_loop;
        r3 = [10*sin(theta_end)+25; 0; h_init - 10*cos(theta_end)];
        r3p = k_loop * [10*cos(theta_end); 0; 10*sin(theta_end)];
        
        % Target for projectile start (adjust these offset values as needed)
        r4 = r3 + [20; 0; 0];
        v4 = [10; 0; 15];
        a0 = r3; a1 = L2 * r3p;
        D = r4 - r3 - a1;
        a3 = L2 * v4 - a1 - 2*D;
        a2 = D - a3;
        r = a0 + a1*s_val2 + a2*s_val2^2 + a3*s_val2^3;
    
    elseif t < T4
        % Segment 5: Projectile (parabolic) motion
        t_par = t - T3p;
        
        % Starting point and velocity from transition target
        theta_end = theta0 + (T3 - T2)*k_loop;
        r4 = [10*sin(theta_end)+35; 0; h_init - 10*cos(theta_end)] + [20; 0; 15];
        v4 = [10; 0; 15];
        r = r4 + v4*t_par + [0; 0; -0.5*g*(t_par.^2)];
    elseif t < T4p
        % Segment 6: Transition from parabolic motion to braking section
        L3 = T4p - T4;
        s_val3 = (t - T4) / L3;
        
        % Get final state from parabolic motion at T4
        t_par_final = T4 - T3p;
        r4_end = r4 + v4 * t_par_final + [0; 0; -0.5 * g * (t_par_final^2)];
        v4_end = v4 + [0; 0; -g * t_par_final];
        
        % Define target state for braking start
        r5 = [r4_end(1) + 100; 0; 0];  % Slightly offset to ensure continuity
        v5 = [max(v4_end(1), 5); 0; 0];  % Ensure a nonzero initial braking velocity
        
        % Cubic transition interpolation
        a0 = r4_end;
        a1 = L3 * v4_end;
        D = r5 - r4_end - a1;
        a3 = L3 * v5 - a1 - 2 * D;
        a2 = D - a3;
        
        r = a0 + a1 * s_val3 + a2 * s_val3^2 + a3 * s_val3^3;
        
    else
        % Segment 7: Braking section (constant deceleration)
        t_brake = t - T4p;
        a_brake = -max(v4_end(1), 5)/t_brake;  % Adjust deceleration as needed
        
        % Initial conditions from transition
        r5 = [r4_end(1) + 100; 0; 0];
        v5 = [max(v4_end(1), 5); 0; 0];
        
        % Apply kinematic equations with constant deceleration
        v_brake = max(0, v5(1) + a_brake * t_brake);
        x_brake = r5(1) + v5(1) * t_brake + 0.5 * a_brake * t_brake^2;
        
        r = [x_brake; 0; 0];

        
    
    end
end

% computeGForces: calculates track arc length, G force, speed, etc.
function [s_track, G_force, x_t, y_t, z_t, v_t, N_dir, Fc_dir] = computeGForces(r_func, h_init, num_points, time_vector, beta_vector)
   
    g = 9.81;
    r_cell = arrayfun(@(tt) r_func(tt), time_vector, 'UniformOutput', false);
    r_t = cell2mat(r_cell);
    x_t = r_t(1,:); y_t = r_t(2,:); z_t = r_t(3,:);
    v_t = sqrt(2 * g * (h_init - z_t));  % energy-based speed
    dt = diff(time_vector);
    dx_dt = diff(x_t)./dt; dy_dt = diff(y_t)./dt; dz_dt = diff(z_t)./dt;
    v_vec = [dx_dt; dy_dt; dz_dt];
    ds = sqrt(dx_dt.^2 + dy_dt.^2 + dz_dt.^2);
    s_track = cumtrapz(time_vector(1:end-1), ds);
    T = v_vec ./ vecnorm(v_vec);
    nT = size(T,2);
    dT_ds = zeros(3, nT);
    
    for i = 2:nT-1
        ds_local = s_track(i+1)-s_track(i-1);
        dT_ds(:,i) = (T(:,i+1)-T(:,i-1)) / ds_local;
    end
   
    dT_ds(:,1) = (T(:,2)-T(:,1)) / (s_track(2)-s_track(1));
    dT_ds(:,end) = (T(:,end)-T(:,end-1)) / (s_track(end)-s_track(end-1));
    curvature = vecnorm(dT_ds);
    n_curv = zeros(size(T));
    
    for i = 1:nT
        if curvature(i) > 1e-6
            n_curv(:,i) = dT_ds(:,i)/curvature(i);
        else
            n_curv(:,i) = [0;0;0];
        end
    end
    
    Fc_dir = n_curv;
    k = [0;0;1];
    N0 = zeros(size(T));
    
    for i = 1:nT
        proj = k - dot(k, T(:,i))*T(:,i);
        
        if norm(proj) < 1e-6
            proj = [1;0;0]; 
        end
        
        N0(:,i) = proj/norm(proj);
    end
    
    N_dir = zeros(size(T));
    beta = deg2rad(beta_vector(1:nT));
    
    for i = 1:nT
        e1 = N0(:,i);
        if norm(n_curv(:,i)) > 1e-6
            e2 = n_curv(:,i) - dot(n_curv(:,i), e1)*e1;
            if norm(e2) > 1e-6
                e2 = e2/norm(e2);
            else
                e2 = -cross(T(:,i), e1); e2 = e2/norm(e2);
            end
        else
            e2 = -cross(T(:,i), e1); e2 = e2/norm(e2);
        end
        N_dir(:,i) = e1*cos(beta(i)) + e2*sin(beta(i));
    end

    v_seg = (v_t(1:end-1) + v_t(2:end)) / 2;
    a_c = v_seg.^2 .* curvature;
    g_vec = [0;0;-g];
    g_comp = zeros(1,nT);

    for i = 1:nT 
        g_comp(i) = -dot(g_vec, N_dir(:,i)); 
    end

    a_N = a_c + g_comp;
    G_force = a_N / g;
    N_dir = [N_dir, N_dir(:,end)];
    Fc_dir = [Fc_dir, Fc_dir(:,end)];
    G_force = [G_force, G_force(end)];
    s_track = [s_track, s_track(end)];
end

% create another function for G-force 
function G_force = gforce (t_values, T1, T2, T3, T3p, T4, T4p, theta0, k_loop, v_t, g)
    % G-force calculation based on different track segments
    % Get the time index based on the time t
    % time_index = find(t_values, 1);  % Find the index where t is located
    
    if t_values < T1 % segment 1: helix loop
        G_force_helix_ud = 1 / cos(theta0); % Example formula for helix
        G_force_helix_fb = 0;
        G_force_helix_lr = 0;

    elseif t_values < T2 % transition segment from helix to loop
        G_force_trans1_ud = 1;
        G_force_trans1_fb = 0;
        G_force_trans1_lr = 0;

    elseif t_values < T3 % loop segment
        theta = theta0 + (t - T2) * k_loop;
        G_force_loop_ud = tan(theta);
        G_force_loop_lr = (v_t(time_index)^2) / (g * 21); % r = 21m
        G_force_loop_fb = 0;

    elseif t_values < T3p % transition segment from loop to parabola
        theta_end = theta0 + (T3 - T2) * k_loop;
        G_force_trans2_ud = tan(theta_end);
        G_force_trans2_lr = (v_t(time_index)^2) / (g * 21);
        G_force_trans2_fb = 0;

    elseif t_values < T4 % parabola segment
        G_force_parabola_ud = 0;
        G_force_parabola_lr = 0;
        G_force_parabola_fb = 0;

    elseif t_values < T4p % transition section from parabola to braking section
        G_force_trans3_ud = 0;
        G_force_trans3_lr = 0;
        G_force_trans3_fb = 0;

    else % braking section segment
        G_force_break_ud = 0;
        G_force_break_lr = 0;
        G_force_break_fb = 0;
    end

    % Store calculated G-forces for each direction
    G_force_fb = [G_force_helix_fb; G_force_trans1_fb; G_force_loop_fb; G_force_trans2_fb; G_force_parabola_fb; G_force_trans3_fb; G_force_break_fb];
    G_force_ud = [G_force_helix_ud; G_force_trans1_ud; G_force_loop_ud; G_force_trans2_ud; G_force_parabola_ud; G_force_trans3_ud; G_force_break_ud];
    G_force_lr = [G_force_helix_lr; G_force_trans1_lr; G_force_loop_lr; G_force_trans2_lr; G_force_parabola_lr; G_force_trans3_lr; G_force_break_lr];

    % Combine all G-forces into a single output
    G_force = [G_force_fb; G_force_ud; G_force_lr];
end