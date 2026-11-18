Mp = 0.127;   % Penduum mass
Lr = 0.216;   % Rotary arm length
Lp = 0.337;   % Pendulum length
Jr = 9.98e-4; % Rotary arm moment of inertia
Jp = 1.2e-2;  % Pendulum moment of inertia
g = 9.81;     % Gravitational acceleration
Kg = 70;      % Total gear ratio
Kt = 7.68e-3; % Motor current-torque constant
Km = 7.68e-3; % Motor back-emf constant
D = 0.0024;   % Viscous damping coefficient
R = 2.6;      % Motor armature resistance

JT = Jp*Mp*Lr^2 + Jr*Jp + 0.25*Jr*Mp*Lp^2;

A = [0, 0, 1, 0;
     0, 0, 0, 1; 
     0, (1/4*Mp^2*Lp^2*Lr*g)/JT, -((Jp + 1/4*Mp*Lp^2)*D)/JT, -(1/2*Mp*Lp*Lr*D)/JT;
     0, (1/2*Mp*Lp*g*(Jr + Mp*Lr^2))/JT, (1/2*Mp*Lp*Lr*D)/JT, -((Jr + Mp*Lr^2)*D)/JT];

B = [0;
     0;
     (Jp + 0.25*Mp*Lp^2)/JT;
     -(0.5*Mp*Lp*Lr)/JT]*Kg*Kt/R;

C = [1, 0, 0, 0;
     0, 1, 0, 0];

D = 0;

del_t = 0.1;   % sampling time

cont_sys = ss(A, B, C, D);   % continuous state space model

discrete_sys = c2d(cont_sys, del_t);  % discrete state space model

I = eye(size(A));
%Ad = expm(A * del_t);
% Bd = (Ad - I)*A^(-1)*B;
Ad = discrete_sys.A;
Bd = discrete_sys.B;
Cd = discrete_sys.C;

% F = [Cd*Ad;
%      Cd*Ad^2;
%      Cd*Ad^3; 
%      Cd*Ad^4; 
%      Cd*Ad^5];
% 
% z = zeros(size(C*Bd));
% 
% Phi = [C*Bd,       z,          z,          z,        z; 
%        C*Ad*Bd,    C*Bd,       z,          z,        z; 
%        C*Ad^2*Bd,  C*Ad*Bd,    C*Bd,       z,        z; 
%        C*Ad^3*Bd,  C*Ad^2*Bd,  C*Ad*Bd,    C*Bd,     z; 
%        C*Ad^4*Bd,  C*Ad^3*Bd,  C*Ad^2*Bd,  C*Ad*Bd,  C*Bd];



% Creating MPC Variables
Np = 15; 
nx = size(Ad, 1);          % Number of states
nu = size(Bd, 2);          % Number of inputs
ny = size(Cd, 1);          % Number of outputs

% Forming the F matrix
s1=size(Cd*Ad);
F=zeros(Np*s1(1),s1(2));
for i=1:Np
    k=(s1(1)*i)-1;
    temp=Cd*Ad^i;
    for j=1:s1(1);
        F(k,:)=temp(j,:);
        k=k+1;
    end
end

s2=size(Cd*Bd);
Phi=zeros(Np*s2(1),Np*s2(2));
for i=1:Np    
    for j=1:i
        temp1=Cd*Ad^(i-j)*Bd;
        Phi(s2(1)*(i-1)+1:s2(1)*i,j)=temp1;
    end
end


u_min = -10;
u_max = 10;
lb = u_min * ones(Np, 1);
ub = u_max * ones(Np, 1);

[m,n]=size(Phi);
q=10;
Qbar = q*eye(m,m);


R=0.1;
Rbar = R*eye(Np);

x0 = [0; 0.5; 0; 0];

%H = Phi' * Qbar * Phi + Rbar;
H = 2 * (Phi' * Qbar * Phi + Rbar);

H=(H+H')/2;

% MPC Simulations
Tsim = 5;
Nsim = Tsim/del_t;

% matrices for logging x(k), y(k) and u(k)
X = zeros(nx, Nsim+1);
Y = zeros(ny, Nsim+1);
U = zeros(nu, Nsim);
uopt=zeros(15,50);

x = x0;

X(:,1) = x;

options = optimoptions('quadprog', 'Display', 'off');

% MPC iteration
for K = 1:Nsim

    %Yfree = F * x;

    %f = 2 * Phi' * Q * (Yfree - R);

    f = 2*Phi'* Qbar*(F*x - R);

    [Uopt, ~, exitflag] = quadprog( ...
        H, f, [], [], [], [], lb, ub, [], options);

    if exitflag <= 0 || isempty(Uopt)
        warning('MPC optimization failed at simulation step %d.', k);
        U(:,k:end) = NaN;
        X(:,k+1:end) = NaN;
        Y(:,k+1:end) = NaN;
        break;
    end
    u = Uopt(1);
    U(:,K) = u;
    uopt(:,K)=Uopt;

    x = Ad*x + Bd*u;

    X(:,K+1) = x;
    Y(:,K+1) = Cd*x;
end

t = 0:del_t:Tsim;


%arm angle 

figure;

plot(t, Y(1,:) * 180/pi, 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Arm Angle (degrees)');
title('Rotary Arm Angle');
grid on;




%pendulum angle


figure;

plot(t, Y(2,:) * 180/pi, 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Pendulum Angle (degrees)');
title('Pendulum Angle');
grid on;





%motor voltages

figure;

stairs(t(1:end-1), U(1,:), 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Motor Voltage (V)');
title('MPC Control Input');
grid on;








