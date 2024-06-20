

rb = 0.22774;   % Fixed platform radius in m
rt = 0.186;   % Moving platform radius in m
gmb = 0.406; % Base angle in radians
gmt = 0.109; % Top angle in radians
z = 0.16;    % Height of the moving platform from base in m
x = 0;    % X-position of the platform in m
     % Y-position of the platform in m
lst = 0.182;  % Length of the struts in m
lc = 0.0125;    % Length of the cranks in m


geom=[rb;rt;gmb;gmt;lst;lc];

q=rss_ikin2(geom);
dq=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

[t,y]=ode45(@diffeq3,[0 0.5],[q;dq]);
plot(y);


rho=2700;% density of crank and strut in kg/m3
den=940;%density of plastic(polymer)in kg/m3
D=0.004;%diameter of crank and strut in m
tk=0.001;%thickness of moving platform in m

mck=(rho*pi*square(D)*lc)/4;
mst=(rho*pi*square(D)*lst)/4;
mp=pi*square(rt)*tk*den;

Ixx1=(1/12)*mck*square(lc);
Iyy1=(1/12)*mck*square(lc);
Izz1=((1/2)*mck*square(D))/4;
Ixx=(1/4)*mp*square(rt);
Iyy=(1/4)*mp*square(rt);
Izz=(1/2)*mp*square(rt)
Ixx2=(1/12)*mst*square(lst);
Iyy2=(1/12)*mst*square(lst);
Izz2=(0.5*mst*square(D))/4;

dim=[rho;den;D;tk;mck;mst;mp];
inertia=[Ixx1;Iyy1;Izz1;Ixx;Iyy;Izz;Ixx2;Iyy2;Izz2];
 TE1=rss_TE(geom,y,dim,inertia);
 ani=duplicate_animation(y);

