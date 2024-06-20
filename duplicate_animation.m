
function ani=duplicate_animation(q)

% % % Platform and Joint Parameters
rb = 22.774;   % Fixed platform radius in cm
rt = 18.6;   % Moving platform radius in cm
gmb = 0.406; % Base angle in radians
gmt = 0.109; % Top angle in radians
z = 16.0;    % Height of the moving platform from base in cm
x = 0;    % X-position of the platform in cm
y = 0;     % Y-position of the platform in cm
lst = 18.2;  % Length of the struts in cm
lc = 1.25;    % Length of the cranks in cm

% % Rotation Parameters
alpha = 0;      % Rotation about x-axis in radians
beta = 0;      % Rotation about y-axis in radians
gamma= 0;      % Rotation about z-axis in radians

geom=[rb;rt;gmb;gmt;z;x;y;lst;lc;alpha;beta;gamma]

%% definition of rotation matrix%%%%%%%%%%%
Rxalpha=[1 0 0;0 cos(alpha) -sin(alpha);0 sin(alpha) cos(alpha)]
Rybeta=[cos(beta) 0 sin(beta);0 1 0;-sin(beta) 0 cos(beta)]
Rzgamma=[cos(gamma) -sin(gamma) 0;sin(gamma) cos(gamma) 0;0 0 1]
Rxyz= Rxalpha*Rybeta*Rzgamma;

%% Angles btw the line (joining each vertices of fixed flatform) and x axis in global frame of reference
gmb1=-gmb;
gmb2=gmb;
gmb3=2*pi/3-gmb;
gmb4=2*pi/3+gmb;
gmb5=4*pi/3-gmb;
gmb6=4*pi/3+gmb;


%% vertices of the fixed platform(bi)w.r.t global frame of reference%%
b=zeros(3,6);

b(:,1)=rb*[cos(-gmb);sin(-gmb);0]
b(:,2)=rb*[cos(gmb);sin(gmb);0]
b(:,3)=rb*[cos(2*pi/3-gmb);sin(2*pi/3-gmb);0]
b(:,4)=rb*[cos(2*pi/3+gmb);sin(2*pi/3+gmb);0]
b(:,5)=rb*[cos(4*pi/3-gmb);sin(4*pi/3-gmb);0]
b(:,6)=rb*[cos(4*pi/3+gmb);sin(4*pi/3+gmb);0]

%% vertices of the moving platform(ti)w.r.t local frame of reference%%
t=zeros(3,6);

t(:,1)=rt*[cos(-gmt);sin(-gmt);0]
t(:,2)=rt*[cos(gmt);sin(gmt);0]
t(:,3)=rt*[cos(2*pi/3-gmt);sin(2*pi/3-gmt);0]
t(:,4)=rt*[cos(2*pi/3+gmt);sin(2*pi/3+gmt);0]
t(:,5)=rt*[cos(4*pi/3-gmt);sin(4*pi/3-gmt);0]
t(:,6)=rt*[cos(4*pi/3+gmt);sin(4*pi/3+gmt);0]



%% tiprime is the new co-ordinate of vertices of the moving platform
%%  w.r.t global frame of reference
tiprime=zeros(3,6)
P=[x;y;z]
for j= 1:6
tiprime(:,j)=P+Rxyz*t(:,j)
end





%% rotation matrix about z %%
    Rz1=[cos(gmb1) -sin(gmb1) 0 ; sin(gmb1) cos(gmb1) 0 ; 0 0 1]
    Rz2=[cos(gmb2) -sin(gmb2) 0;sin(gmb2) cos(gmb2) 0;0 0 1]
    Rz3=[cos(gmb3) -sin(gmb3) 0;sin(gmb3) cos(gmb3) 0;0 0 1]
    Rz4=[cos(gmb4) -sin(gmb4) 0;sin(gmb4) cos(gmb4) 0;0 0 1]
    Rz5=[cos(gmb5) -sin(gmb5) 0;sin(gmb5) cos(gmb5) 0;0 0 1]
    Rz6=[cos(gmb6) -sin(gmb6) 0;sin(gmb6) cos(gmb6) 0;0 0 1]




        th1= q(1);
        th2= q(2);
        th3= q(3);
        th4= q(4);
        th5= q(5);
        th6= q(6);
% Calculate vertices of the cranks (ck) using rotation matrices
ck = zeros(3, 6);
ck(:,1) = Rz1 * [rb; lc * cos(th1); lc * sin(th1)];
ck(:,2) = Rz2 * [rb; lc * cos(th2); lc * sin(th2)];
ck(:,3) = Rz3 * [rb; lc * cos(th3); lc * sin(th3)];
ck(:,4) = Rz4 * [rb; lc * cos(th4); lc * sin(th4)];
ck(:,5) = Rz5 * [rb; lc * cos(th5); lc * sin(th5)];
ck(:,6) = Rz6 * [rb; lc * cos(th6); lc * sin(th6)];


% Plotting the Stewart Platform
figure;
hold on;

% Plot fixed platform (bi)
plot3([b(1,1), b(1,2), b(1,3), b(1,4), b(1,5), b(1,6), b(1,1)], ...
      [b(2,1), b(2,2), b(2,3), b(2,4), b(2,5), b(2,6), b(2,1)], ...
      [b(3,1), b(3,2), b(3,3), b(3,4), b(3,5), b(3,6), b(3,1)], 'b-', 'LineWidth', 2); % Fixed platform

% Plot moving platform (tiprime)
plot3([tiprime(1,1), tiprime(1,2), tiprime(1,3), tiprime(1,4), tiprime(1,5), tiprime(1,6), tiprime(1,1)], ...
      [tiprime(2,1), tiprime(2,2), tiprime(2,3), tiprime(2,4), tiprime(2,5), tiprime(2,6), tiprime(2,1)], ...
      [tiprime(3,1), tiprime(3,2), tiprime(3,3), tiprime(3,4), tiprime(3,5), tiprime(3,6), tiprime(3,1)], 'r-', 'LineWidth', 2); % Moving platform

% Plot cranks (ck) connecting to fixed platform vertices (b)
for j = 1:6
    plot3([b(1,j), ck(1,j)], [b(2,j), ck(2,j)], [b(3,j), ck(3,j)], 'k-', 'LineWidth', 1.5); % Connect to fixed platform
    plot3([ck(1,j), tiprime(1,j)], [ck(2,j), tiprime(2,j)], [ck(3,j), tiprime(3,j)], 'g-', 'LineWidth', 1.5); % Connect to moving platform
end

% Additional settings for visualization
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
grid on;
axis equal; % Ensure equal scaling for all axes
view(3); % Set the 3D view
title('Stewart Platform Visualization');

hold off;
