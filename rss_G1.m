function G=rss_G1(geom,inertia,q,dim);
Ixx1=inertia(1);
Iyy1=inertia(2);
Izz1=inertia(3);
Ixx=inertia(4);
Iyy=inertia(5);
Izz=inertia(6);
Ixx2=inertia(7);
Iyy2=inertia(8);
Izz2=inertia(9);
mck=dim(5);
mst=dim(6);
mp=dim(7);

rb = geom(1);
rt = geom(2);
gmb= geom(3);
gmt = geom(4);
lst = geom(5);
lc = geom(6);

%q=[0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5]
%dq=[0.1;0;0;0.1;0.5;0;0.5;0.5;0.2;0.5;0.5;0.5;0.5;0.5;0.5;0;0.5;0.5;0.5;0.5;0.5;0;0.5;0]
%q=[0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5]

th1= q(1);
th2= q(2);
th3= q(3);
th4= q(4);
th5= q(5);
th6= q(6);

phi1= q(7);
phi2= q(8);
phi3= q(9);
phi4= q(10);
phi5= q(11);
phi6= q(12);

shi1= q(13);
shi2= q(14);
shi3= q(15);
shi4= q(16);
shi5= q(17);
shi6= q(18);

x= q(19);
y= q(20);
z= q(21);
a= q(22);
b= q(23);
c= q(24);

G=[(49*lc*mck*cos(th1))/10 + (49*lc*mst*cos(th1))/5, (49*lc*mck*cos(th2))/10 + (49*lc*mst*cos(th2))/5, (49*lc*mck*cos(th3))/10 + (49*lc*mst*cos(th3))/5, (49*lc*mck*cos(th4))/10 + (49*lc*mst*cos(th4))/5, (49*lc*mck*cos(th5))/10 + (49*lc*mst*cos(th5))/5, (49*lc*mck*cos(th6))/10 + (49*lc*mst*cos(th6))/5, -(49*lst*mst*cos(shi1)*sin(phi1))/10, -(49*lst*mst*cos(shi2)*sin(phi2))/10, -(49*lst*mst*cos(shi3)*sin(phi3))/10, -(49*lst*mst*cos(shi4)*sin(phi4))/10, -(49*lst*mst*cos(shi5)*sin(phi5))/10, -(49*lst*mst*cos(shi6)*sin(phi6))/10, -(49*lst*mst*cos(phi1)*sin(shi1))/10, -(49*lst*mst*cos(phi2)*sin(shi2))/10, -(49*lst*mst*cos(phi3)*sin(shi3))/10, -(49*lst*mst*cos(phi4)*sin(shi4))/10, -(49*lst*mst*cos(phi5)*sin(shi5))/10, -(49*lst*mst*cos(phi6)*sin(shi6))/10, 0, 0, (49*mp)/5, 0, 0, 0];
