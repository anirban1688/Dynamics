%%%%%% Symbolic derivation of the equation of motion of the spring pendulum

clc
syms l r(t) theta(t) m g t rdot thtdot ath ar k SS Y  %%%% Define the symbolic variables

x=(l+r)*sin(theta); %%%% get the position and velocity of the bob
y=(l+r)*cos(theta);
vx=diff(x,t);
vy=diff(y,t);
vx1=subs(vx,[diff(r(t), t),diff(theta(t), t)],[rdot,thtdot]);
vy1=subs(vy,[diff(r(t), t),diff(theta(t), t)],[rdot,thtdot]);

T=simplify((1/2)*m*(vx1^2+vy1^2)); %%%%%% Calculate the kinetic and potential energy
V=m*g*y + (1/2)*k*r^2;

L=T-V; %%%%%% Calculate the Lagrangian 

%%%%%%% Derivation of the equation of motion
%%%%%%% computation of d/dt(del L/del thetadot);
temp1a=diff(L,thtdot);
temp1b=subs(temp1a,[rdot,thtdot],[diff(r(t), t),diff(theta(t), t)]);
temp1c=diff(temp1b,t);
temp1d=subs(temp1c,[diff(theta(t), t, t),diff(r(t), t, t)],[ath,ar]);
eq1a=simplify(subs(temp1d,[diff(r(t), t),diff(theta(t), t)],[rdot,thtdot]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eq1b=diff(L,theta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eqn1=simplify((eq1a-eq1b)/(m*(l+r)))
%%eqn1a=simplify(subs(eqn1,[ar,ath,rdot,thtdot],[diff(r(t), t, t),diff(theta(t), t, t),diff(r(t), t),diff(theta(t), t)]));

%%%%%%%% computation of d/dt(del L/del rdot)%%%%%%%%%%%%%%%%%%%%%
temp2a=diff(L,rdot);
temp2b=subs(temp2a,[rdot,thtdot],[diff(r(t), t),diff(theta(t), t)]);
temp2c=diff(temp2b,t);
temp2d=subs(temp2c,[diff(theta(t), t, t),diff(r(t), t, t)],[ath,ar]);
eq2a=simplify(subs(temp2d,[diff(r(t), t),diff(theta(t), t)],[rdot,thtdot]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eq2b=diff(L,r);
eqn2=eq2a-eq2b
%%eqn2a=simplify(subs(eqn2,[ar,ath,rdot,thtdot],[diff(r(t), t, t),diff(theta(t), t, t),diff(r(t), t),diff(theta(t), t)]));

%%%%% Solve for the accelerations %%%%%%%%%%
temp3=solve(eqn1==0,ath)
temp4=solve(eqn2==0,ar)

%%%%% Converting to state-space and solving the differential equations
%%%%% using ode45

Eq1 = diff(theta,2) == subs(temp3,[thtdot,rdot],[diff(theta),diff(r)])
Eq2 = diff(r,2) == subs(temp4,[thtdot,rdot],[diff(theta),diff(r)])


[VF,Sbs] = odeToVectorField(Eq1,Eq2)
spring_pendulum = matlabFunction(VF, 'Vars',{t,Y,g,l,m,k})

g=9.8;
k=0.5;
m=1;
l=1;

tspan = [0,20];t0 = [0.1,0,pi/6,0];
[t,v_z] = ode45(@(T,Y)spring_pendulum(t,Y,g,l,m,k),tspan,t0);
figure
plot(t,v_z(:,1))
grid

title('Solution of simple pendulum with ODE45');
xlabel('Time t');
ylabel('State variables');
legend('r','$\dot{r}$','$\theta$','$\dot{\theta}$','Interpreter','latex')


