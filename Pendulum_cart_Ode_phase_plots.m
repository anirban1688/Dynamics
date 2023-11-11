%%%%%% Symbolic derivation of the equation of motion of the pendulum on cart
%NOTE1: I'm taking "z(t)" instead of "x(t)" to avoid any discrepency between "x(t)" as variable and x-axis components.
%NOTE2: az = (z)doubledot and ath = (theta)doubledot
clc
clear all
syms l b M z(t) theta(t) m g t zdot thtdot ath az Y u %%%% Define the symbolic variables 

x=z+l*sin(theta); %%%% get the position and velocity of the bob
y=b+l*cos(theta);
vx=diff(x,t);
vy=diff(y,t);
vx1=subs(vx,[diff(z(t), t),diff(theta(t), t)],[zdot,thtdot]);
vy1=subs(vy,[diff(z(t), t),diff(theta(t), t)],[zdot,thtdot]);

T=simplify((1/2)*m*(vx1^2+vy1^2)+(1/2)*M*(zdot)^2); %%%%%% Calculate the kinetic and potential energy
V=m*g*(b/2) + m*g*(b+l*cos(theta));

L=T-V; %%%%%% Calculate the Lagrangian 

%%%%%%% Derivation of the equation of motion
%%%%%%% computation of d/dt(del L/del thetadot);
temp1a=diff(L,thtdot);
temp1b=subs(temp1a,[zdot,thtdot],[diff(z(t), t),diff(theta(t), t)]);
temp1c=diff(temp1b,t);
temp1d=subs(temp1c,[diff(theta(t), t, t),diff(z(t), t, t)],[ath,az]);
eq1a=simplify(subs(temp1d,[diff(z(t), t),diff(theta(t), t)],[zdot,thtdot]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eq1b=diff(L,theta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eqn1=simplify((eq1a-eq1b)/m)
%%eqn1a=simplify(subs(eqn1,[ar,ath,rdot,thtdot],[diff(r(t), t, t),diff(theta(t), t, t),diff(r(t), t),diff(theta(t), t)]));

%%%%%%%% computation of d/dt(del L/del rdot)%%%%%%%%%%%%%%%%%%%%%
temp2a=diff(L,zdot);
temp2b=subs(temp2a,[zdot,thtdot],[diff(z(t), t),diff(theta(t), t)]);
temp2c=diff(temp2b,t);
temp2d=subs(temp2c,[diff(theta(t), t, t),diff(z(t), t, t)],[ath,az]);
eq2a=simplify(subs(temp2d,[diff(z(t), t),diff(theta(t), t)],[zdot,thtdot]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eq2b=diff(L,z);
eqn2=eq2a-eq2b
%%eqn2a=simplify(subs(eqn2,[ar,ath,zdot,thtdot],[diff(z(t), t, t),diff(theta(t), t, t),diff(z(t), t),diff(theta(t), t)]));

[temp3,temp4] = solve(eqn1,eqn2 == 0, ath, az)
%%%%% Solve for the accelerations %%%%%%%%%%
% temp3=solve(eqn1==0,ath)
% temp4=solve(eqn2==0,az)
% 
% %%%%% Converting to state-space and solving the differential equations
% %%%%% using ode45
% 
Eq1 = diff(theta,2) == subs(temp3,[thtdot,zdot],[diff(theta),diff(z)])
Eq2 = diff(z,2) == subs(temp4,[thtdot,zdot],[diff(theta),diff(z)])
% 
% 
[VF,Sbs] = odeToVectorField(Eq1,Eq2)
Pendulum_on_cart = matlabFunction(VF, 'Vars',{t,Y,g,l,m,M,b})

g=9.8;
b=1;
M=1;
m=1;
l=1;

tspan = [0,20];t0 = [0.1,0,pi/6,0];
[t,v_z] = ode45(@(T,Y)Pendulum_on_cart(t,Y,g,l,m,M,b),tspan,t0);

q=v_z(:,1)
e=v_z(:,2)
d=v_z(:,3)
c=v_z(:,4)

 f1=figure('Name','Pendulum_on_cart plot','NumberTitle','off');
 plot(t,v_z)
 grid
 
 title('Solution of Pendulum on cart system');
 xlabel('Time t');
 ylabel('State variables');
 legend('z','$\dot{z}$','$\theta$','$\dot{\theta}$','Interpreter','latex')
% 


y1 = linspace(-2,6,30);
y2 = linspace(-6,6,30);


% creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 columns
[x1,y1] = meshgrid(y1,y2);
size(x1)
size(y1)

u1 = zeros(size(x1));
v1 = zeros(size(x1));
u2 = zeros(size(x1));
v2 = zeros(size(x1));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(x1)
    Yprime = Pendulum_on_cart(t,[x1(i); y1(i); x1(i); y1(i)],g,l,m,M,b);
    u1(i) = Yprime(1);
    v1(i) = Yprime(2);
    u2(i) = Yprime(3);
    v2(i) = Yprime(4);
end


figure;
quiver(x1,y1,u1,v1,'r'); %figure(gcf)
figure;
quiver(x1,y1,u2,v2,'b');
hold on
%for y20=[0 0.5 1 1.5 2 2.5]
    [ts,ys] = ode45(@(T,Y)Pendulum_on_cart(t,Y,g,l,m,M,b),[0,10],[0,0,pi/6,0]);
    %[ts,ys]=ode45(@(T,Y)Pendulum_on_cart(t,Y,g,l,m,M,b),tspan1,t0);
    plot(ys(:,3),ys(:,4))
    plot(ys(1,1),ys(1,2),'bo') % starting point
    plot(ys(end,1),ys(end,2),'ks') % ending point
%end
hold off

%%%%%%% Derivation of the forced system and control laws %%%%%%%%%%%%

[temp31,temp41] = solve(eqn1,eqn2-u == 0, ath, az)

Eq11 = diff(theta,2) == subs(temp31,[thtdot,zdot],[diff(theta),diff(z)])
Eq21 = diff(z,2) == subs(temp41,[thtdot,zdot],[diff(theta),diff(z)])

[VF1,Sbs1] = odeToVectorField(Eq11,Eq21)
Pendulum_on_cart_inp = matlabFunction(VF1, 'Vars',{t,Y,g,l,m,M,b,u})