function mass=rss_mass2(geom,inertia,q,dim)



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

mass=[Izz1 + (lc^2*mck)/4 + lc^2*mst, 0, 0, 0, 0, 0, - (lc*lst*mst*cos(shi1 )*cos(th1 )*sin(phi1 ))/2 - (lc*lst*mst*cos(phi1 )*cos(shi1 )*sin(th1 )*sin(gmb))/2, 0, 0, 0, 0, 0, (lc*lst*mst*cos(shi1 )*sin(th1 )*cos(gmb))/2 - (lc*lst*mst*cos(phi1 )*cos(th1 )*sin(shi1 ))/2 + (lc*lst*mst*sin(phi1 )*sin(shi1 )*sin(th1 )*sin(gmb))/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, Izz1 + (lc^2*mck)/4 + lc^2*mst, 0, 0, 0, 0, 0, (lc*lst*mst*cos(phi2 )*cos(shi2 )*sin(th2 )*sin(gmb))/2 - (lc*lst*mst*cos(shi2 )*cos(th2 )*sin(phi2 ))/2, 0, 0, 0, 0, 0, (lc*lst*mst*cos(shi2 )*sin(th2 )*cos(gmb))/2 - (lc*lst*mst*cos(phi2 )*cos(th2 )*sin(shi2 ))/2 - (lc*lst*mst*sin(phi2 )*sin(shi2 )*sin(th2 )*sin(gmb))/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, Izz1 + (lc^2*mck)/4 + lc^2*mst, 0, 0, 0, 0, 0, (lc*lst*mst*cos(phi3 )*cos(shi3 )*sin(th3 )*sin(gmb))/4 - (lc*lst*mst*cos(shi3 )*cos(th3 )*sin(phi3 ))/2 + (3^(1/2)*lc*lst*mst*cos(phi3 )*cos(shi3 )*sin(th3 )*cos(gmb))/4, 0, 0, 0, 0, 0, (3^(1/2)*lc*lst*mst*cos(shi3 )*sin(th3 )*sin(gmb))/4 - (lc*lst*mst*cos(phi3 )*cos(th3 )*sin(shi3 ))/2 - (lc*lst*mst*cos(shi3 )*sin(th3 )*cos(gmb))/4 - (lc*lst*mst*sin(phi3 )*sin(shi3 )*sin(th3 )*sin(gmb))/4 - (3^(1/2)*lc*lst*mst*sin(phi3 )*sin(shi3 )*sin(th3 )*cos(gmb))/4, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, Izz1 + (lc^2*mck)/4 + lc^2*mst, 0, 0, 0, 0, 0, (3^(1/2)*lc*lst*mst*cos(phi4 )*cos(shi4 )*sin(th4 )*cos(gmb))/4 - (lc*lst*mst*cos(phi4 )*cos(shi4 )*sin(th4 )*sin(gmb))/4 - (lc*lst*mst*cos(shi4 )*cos(th4 )*sin(phi4 ))/2, 0, 0, 0, 0, 0, (lc*lst*mst*sin(phi4 )*sin(shi4 )*sin(th4 )*sin(gmb))/4 - (lc*lst*mst*cos(phi4 )*cos(th4 )*sin(shi4 ))/2 - (3^(1/2)*lc*lst*mst*cos(shi4 )*sin(th4 )*sin(gmb))/4 - (lc*lst*mst*cos(shi4 )*sin(th4 )*cos(gmb))/4 - (3^(1/2)*lc*lst*mst*sin(phi4 )*sin(shi4 )*sin(th4 )*cos(gmb))/4, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, Izz1 + (lc^2*mck)/4 + lc^2*mst, 0, 0, 0, 0, 0, (lc*lst*mst*cos(phi5 )*cos(shi5 )*sin(th5 )*sin(gmb))/4 - (lc*lst*mst*cos(shi5 )*cos(th5 )*sin(phi5 ))/2 - (3^(1/2)*lc*lst*mst*cos(phi5 )*cos(shi5 )*sin(th5 )*cos(gmb))/4, 0, 0, 0, 0, 0, (3^(1/2)*lc*lst*mst*sin(phi5 )*sin(shi5 )*sin(th5 )*cos(gmb))/4 - (lc*lst*mst*cos(phi5 )*cos(th5 )*sin(shi5 ))/2 - (3^(1/2)*lc*lst*mst*cos(shi5 )*sin(th5 )*sin(gmb))/4 - (lc*lst*mst*sin(phi5 )*sin(shi5 )*sin(th5 )*sin(gmb))/4 - (lc*lst*mst*cos(shi5 )*sin(th5 )*cos(gmb))/4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, Izz1 + (lc^2*mck)/4 + lc^2*mst, 0, 0, 0, 0, 0, - (lc*lst*mst*cos(shi6 )*cos(th6 )*sin(phi6 ))/2 - (lc*lst*mst*cos(phi6 )*cos(shi6 )*sin(th6 )*sin(gmb))/4 - (3^(1/2)*lc*lst*mst*cos(phi6 )*cos(shi6 )*sin(th6 )*cos(gmb))/4, 0, 0, 0, 0, 0, (3^(1/2)*lc*lst*mst*cos(shi6 )*sin(th6 )*sin(gmb))/4 - (lc*lst*mst*cos(phi6 )*cos(th6 )*sin(shi6 ))/2 - (lc*lst*mst*cos(shi6 )*sin(th6 )*cos(gmb))/4 + (lc*lst*mst*sin(phi6 )*sin(shi6 )*sin(th6 )*sin(gmb))/4 + (3^(1/2)*lc*lst*mst*sin(phi6 )*sin(shi6 )*sin(th6 )*cos(gmb))/4, 0, 0, 0, 0, 0, 0; - (lc*lst*mst*cos(shi1 )*cos(th1 )*sin(phi1 ))/2 - (lc*lst*mst*cos(phi1 )*cos(shi1 )*sin(th1 )*sin(gmb))/2, 0, 0, 0, 0, 0, Izz2 + Iyy2*cos(shi1 )^2 - Izz2*cos(shi1 )^2 + (lst^2*mst*cos(shi1 )^2)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, (lc*lst*mst*cos(phi2 )*cos(shi2 )*sin(th2 )*sin(gmb))/2 - (lc*lst*mst*cos(shi2 )*cos(th2 )*sin(phi2 ))/2, 0, 0, 0, 0, 0, Izz2 + Iyy2*cos(shi2 )^2 - Izz2*cos(shi2 )^2 + (lst^2*mst*cos(shi2 )^2)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, (lc*lst*mst*cos(phi3 )*cos(shi3 )*sin(th3 )*sin(gmb))/4 - (lc*lst*mst*cos(shi3 )*cos(th3 )*sin(phi3 ))/2 + (3^(1/2)*lc*lst*mst*cos(phi3 )*cos(shi3 )*sin(th3 )*cos(gmb))/4, 0, 0, 0, 0, 0, Izz2 + Iyy2*cos(shi3 )^2 - Izz2*cos(shi3 )^2 + (lst^2*mst*cos(shi3 )^2)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, (3^(1/2)*lc*lst*mst*cos(phi4 )*cos(shi4 )*sin(th4 )*cos(gmb))/4 - (lc*lst*mst*cos(phi4 )*cos(shi4 )*sin(th4 )*sin(gmb))/4 - (lc*lst*mst*cos(shi4 )*cos(th4 )*sin(phi4 ))/2, 0, 0, 0, 0, 0, Izz2 + Iyy2*cos(shi4 )^2 - Izz2*cos(shi4 )^2 + (lst^2*mst*cos(shi4 )^2)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, (lc*lst*mst*cos(phi5 )*cos(shi5 )*sin(th5 )*sin(gmb))/4 - (lc*lst*mst*cos(shi5 )*cos(th5 )*sin(phi5 ))/2 - (3^(1/2)*lc*lst*mst*cos(phi5 )*cos(shi5 )*sin(th5 )*cos(gmb))/4, 0, 0, 0, 0, 0, Izz2 + Iyy2*cos(shi5 )^2 - Izz2*cos(shi5 )^2 + (lst^2*mst*cos(shi5 )^2)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, - (lc*lst*mst*cos(shi6 )*cos(th6 )*sin(phi6 ))/2 - (lc*lst*mst*cos(phi6 )*cos(shi6 )*sin(th6 )*sin(gmb))/4 - (3^(1/2)*lc*lst*mst*cos(phi6 )*cos(shi6 )*sin(th6 )*cos(gmb))/4, 0, 0, 0, 0, 0, Izz2 + Iyy2*cos(shi6 )^2 - Izz2*cos(shi6 )^2 + (lst^2*mst*cos(shi6 )^2)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; (lc*lst*mst*cos(shi1 )*sin(th1 )*cos(gmb))/2 - (lc*lst*mst*cos(phi1 )*cos(th1 )*sin(shi1 ))/2 + (lc*lst*mst*sin(phi1 )*sin(shi1 )*sin(th1 )*sin(gmb))/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Ixx2 + (lst^2*mst)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, (lc*lst*mst*cos(shi2 )*sin(th2 )*cos(gmb))/2 - (lc*lst*mst*cos(phi2 )*cos(th2 )*sin(shi2 ))/2 - (lc*lst*mst*sin(phi2 )*sin(shi2 )*sin(th2 )*sin(gmb))/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Ixx2 + (lst^2*mst)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, (3^(1/2)*lc*lst*mst*cos(shi3 )*sin(th3 )*sin(gmb))/4 - (lc*lst*mst*cos(phi3 )*cos(th3 )*sin(shi3 ))/2 - (lc*lst*mst*cos(shi3 )*sin(th3 )*cos(gmb))/4 - (lc*lst*mst*sin(phi3 )*sin(shi3 )*sin(th3 )*sin(gmb))/4 - (3^(1/2)*lc*lst*mst*sin(phi3 )*sin(shi3 )*sin(th3 )*cos(gmb))/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Ixx2 + (lst^2*mst)/4, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, (lc*lst*mst*sin(phi4 )*sin(shi4 )*sin(th4 )*sin(gmb))/4 - (lc*lst*mst*cos(phi4 )*cos(th4 )*sin(shi4 ))/2 - (3^(1/2)*lc*lst*mst*cos(shi4 )*sin(th4 )*sin(gmb))/4 - (lc*lst*mst*cos(shi4 )*sin(th4 )*cos(gmb))/4 - (3^(1/2)*lc*lst*mst*sin(phi4 )*sin(shi4 )*sin(th4 )*cos(gmb))/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Ixx2 + (lst^2*mst)/4, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, (3^(1/2)*lc*lst*mst*sin(phi5 )*sin(shi5 )*sin(th5 )*cos(gmb))/4 - (lc*lst*mst*cos(phi5 )*cos(th5 )*sin(shi5 ))/2 - (3^(1/2)*lc*lst*mst*cos(shi5 )*sin(th5 )*sin(gmb))/4 - (lc*lst*mst*sin(phi5 )*sin(shi5 )*sin(th5 )*sin(gmb))/4 - (lc*lst*mst*cos(shi5 )*sin(th5 )*cos(gmb))/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Ixx2 + (lst^2*mst)/4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, (3^(1/2)*lc*lst*mst*cos(shi6 )*sin(th6 )*sin(gmb))/4 - (lc*lst*mst*cos(phi6 )*cos(th6 )*sin(shi6 ))/2 - (lc*lst*mst*cos(shi6 )*sin(th6 )*cos(gmb))/4 + (lc*lst*mst*sin(phi6 )*sin(shi6 )*sin(th6 )*sin(gmb))/4 + (3^(1/2)*lc*lst*mst*sin(phi6 )*sin(shi6 )*sin(th6 )*cos(gmb))/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Ixx2 + (lst^2*mst)/4, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, mp, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, mp, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, mp, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Izz + Iyy*cos(b )^2 - Izz*cos(b )^2 + Ixx*cos(b )^2*cos(c )^2 - Iyy*cos(b )^2*cos(c )^2, Ixx*cos(b )*cos(c )*sin(c ) - Iyy*cos(b )*cos(c )*sin(c ), Izz*sin(b ); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Ixx*cos(b )*cos(c )*sin(c ) - Iyy*cos(b )*cos(c )*sin(c ), Ixx - Ixx*cos(c )^2 + Iyy*cos(c )^2, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Izz*sin(b ), 0, Izz];
