clear;
clc;
close all;

params.Mass=1;
params.Lenght=1;
params.g=9.81;
params.b=-1; % damping
params.I=params.Mass*params.Lenght^2;

params.Slope=5;
params.K=10;

x0=[1;0]; % defining theta and omega [theta;omega]

tspan=[0 10]; % Simulation time

[t,x]=ode45(@(t,x) pendulum_smc(t,x,params),tspan,x0);

theta=x(:,1);
omega=x(:,2);

S=zeros(length(t),1);
u=zeros(length(t),1);

for i=1:length(t)

    x1=theta(i);
    x2=omega(i);

    M=params.Mass;
    L=params.Lenght;
    g=params.g;
    b=params.b;
    I=params.I;
    Slope=params.Slope;
    K=params.K;

    % Sliding surface
    S(i)=Slope*x1+x2;

    % Equivalent control
    u_eq=b*x2+M*g*L*sin(x1)-I*Slope*x2;

    % Switching control
    u_sw=-I*K*sign(S(i));

    % Total control
    u(i)=u_eq+u_sw;

end

figure

subplot(2,2,1)
plot(t,theta,'LineWidth',2)
xlabel('Time (s)')
ylabel('\theta (rad)')
title('Pendulum Angle vs Time')
grid on

subplot(2,2,2)
plot(theta,omega,'LineWidth',2)
xlabel('\theta (rad)')
ylabel('\omega (rad/s)')
title('Phase Portrait')
grid on

subplot(2,2,3)
plot(t,u,'LineWidth',2)
xlabel('Time (s)')
ylabel('Control Input (Torque)')
title('Control Input')
grid on

subplot(2,2,4)
plot(t,S,'LineWidth',2)
xlabel('Time (s)')
ylabel('Sliding Surface S')
title('Sliding Surface Convergence')
grid on


function dxdt=pendulum_smc(~,x,params)

x1=x(1); % theta
x2=x(2); % theta_dot

M=params.Mass;
L=params.Lenght;
g=params.g;
b=params.b;
I=params.I;
Slope=params.Slope;
K=params.K;

% Sliding surface
S=Slope*x1+x2;

% Equivalent control
u_eq=b*x2+M*g*L*sin(x1)-I*Slope*x2;

% Switching control
u_sw=-I*K*sign(S);

% Total control
u=u_eq+u_sw;

% System dynamics
dx1=x2;
dx2=-(b/I)*x2-(M*g*L/I)*sin(x1)+(1/I)*u;

dxdt=[dx1;dx2];

end