clc;
clear;
close all;

% Parameters
m1 = 1; l1 = 1; me = 2; delta_e = deg2rad(30); % converting 30° to radians
I1 = 0.12; lc1 = 0.5; le = 0.25; lce = 0.6;

% Constants derived from parameters
a1 = I1 + m1 * lc1^2 + le + me * lce^2 + me * l1^2;
a2 = le + me * lce^2;
a3 = me * l1 * lce * cos(delta_e);
a4 = me * l1 * lce * sin(delta_e);

% Initial conditions with disturbances
q1 = 0 ; 
q2 = 0; % adding disturbance to initial positions
q1_dot = 0 ; 
q2_dot = 0 ; % adding disturbance to initial velocities
q1_desired = deg2rad(60); % desired position q1 in radians
q2_desired = deg2rad(90); % desired position q2 in radians

% Time parameters
tspan = [0 5]; % simulation duration

% Control gains
Kd = 100; % derivative gain
Kp = 20 * Kd; % proportional gain

% Solve using ode45
[t, y] = ode45(@(t, y) dynamics(t, y, t, Kp, Kd, q1_desired, q2_desired, a1, a2, a3, a4), tspan, [q1; q2; q1_dot; q2_dot]);

% Extract results
q1_array = y(:, 1);
q2_array = y(:, 2);

% Compute position errors
eq1_error = rad2deg(q1_array - q1_desired);
eq2_error = rad2deg(q2_array - q2_desired);

% Compute control torques directly during ODE integration
tau_history = zeros(length(t), 2);
for i = 1:length(t)
    [~, tau1, tau2] = dynamics(t(i), y(i, :)', t(i), Kp, Kd, q1_desired, q2_desired, a1, a2, a3, a4);
    tau_history(i, :) = [tau1, tau2];
end

% Plot results
figure;
subplot(2,1,1);
plot(t, eq1_error, 'r', 'LineWidth', 1.5);
hold on;
plot(t, eq2_error, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Position Errors (deg)');
legend('q1 error', 'q2 error');
grid on;

subplot(2,1,2);
plot(t, tau_history(:,1), 'r', 'LineWidth', 1.5);
hold on;
plot(t, tau_history(:,2), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Torques (Nm)');
legend('\tau_1', '\tau_2');
grid on;

% Function defining system dynamics
function [dydt, tau1, tau2] = dynamics(~, y, t, Kp, Kd, q1_desired, q2_desired, a1, a2, a3, a4)

    u_lim=pi/100;
    l_lim=pi/200;
    di=l_lim+(u_lim-l_lim)*rand(1,1);
    q1 = y(1)+di;
    q2 = y(2)-di;
    q1_dot = y(3);
    q2_dot = y(4);

    H = [a1 + 2 * a3 * cos(q2) + 2 * a4 * sin(q2), a2 + a3 * cos(q2) + a4 * sin(q2);
         a2 + a3 * cos(q2) + a4 * sin(q2),        a2];
    
    h = a3 * sin(q2) - a4 * cos(q2);

    u_lim=pi/100;
    l_lim=pi/200;
    di=l_lim+(u_lim-l_lim)*rand(1,1);
    
    q1_err = q1_desired - q1;
    q2_err = q2_desired - q2;
    q1_dot_err = -q1_dot;
    q2_dot_err = -q2_dot;
    
    
    % Time-varying disturbance
    %disturbance = [5 * sin(2 * pi * t); 5 * cos(2 * pi * t)*sin(2*pi*t)];
    
    tau=[Kp * q1_err + Kd * q1_dot_err;
         Kp * q2_err + Kd * q2_dot_err];
    
    tau1 = tau(1);
    tau2 = tau(2);
    
    q_ddot = H \ tau;  % Solve for q1_ddot and q2_ddot

    dydt = [q1_dot; q2_dot; q_ddot];
end
