function J=rss_jac2(geom,q,dim)



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

J=[lc*sin(th1)*sin(gmb), 0, 0, 0, 0, 0, -lst*cos(phi1)*cos(shi1), 0, 0, 0, 0, 0, lst*sin(phi1)*sin(shi1), 0, 0, 0, 0, 0, 0, 0, 0, 0, - rt*cos(c)*sin(b)*cos(gmt) - rt*sin(b)*sin(c)*sin(gmt), rt*cos(b)*cos(c)*sin(gmt) - rt*cos(b)*sin(c)*cos(gmt); lc*sin(th1)*cos(gmb), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, lst*cos(shi1), 0, 0, 0, 0, 0, 0, 0, 0, rt*sin(gmt)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - rt*cos(gmt)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)), rt*cos(b)*cos(c)*sin(a)*cos(gmt) + rt*cos(b)*sin(a)*sin(c)*sin(gmt), rt*cos(gmt)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) + rt*sin(gmt)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)); -lc*cos(th1), 0, 0, 0, 0, 0, lst*cos(shi1)*sin(phi1), 0, 0, 0, 0, 0, lst*cos(phi1)*sin(shi1), 0, 0, 0, 0, 0, 0, 0, 0, rt*cos(gmt)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) - rt*sin(gmt)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)), - rt*cos(a)*cos(b)*cos(c)*cos(gmt) - rt*cos(a)*cos(b)*sin(c)*sin(gmt), rt*cos(gmt)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) + rt*sin(gmt)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)); 0, -lc*sin(th2)*sin(gmb), 0, 0, 0, 0, 0, -lst*cos(phi2)*cos(shi2), 0, 0, 0, 0, 0, lst*sin(phi2)*sin(shi2), 0, 0, 0, 0, 0, 0, 0, 0, rt*sin(b)*sin(c)*sin(gmt) - rt*cos(c)*sin(b)*cos(gmt), - rt*cos(b)*cos(c)*sin(gmt) - rt*cos(b)*sin(c)*cos(gmt); 0, lc*sin(th2)*cos(gmb), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, lst*cos(shi2), 0, 0, 0, 0, 0, 0, 0, - rt*cos(gmt)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) - rt*sin(gmt)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)), rt*cos(b)*cos(c)*sin(a)*cos(gmt) - rt*cos(b)*sin(a)*sin(c)*sin(gmt), rt*cos(gmt)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - rt*sin(gmt)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)); 0, -lc*cos(th2), 0, 0, 0, 0, 0, lst*cos(shi2)*sin(phi2), 0, 0, 0, 0, 0, lst*cos(phi2)*sin(shi2), 0, 0, 0, 0, 0, 0, 0, rt*cos(gmt)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) + rt*sin(gmt)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)), rt*cos(a)*cos(b)*sin(c)*sin(gmt) - rt*cos(a)*cos(b)*cos(c)*cos(gmt), rt*cos(gmt)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - rt*sin(gmt)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)); 0, 0, lc*sin(th3)*sin(gmb - (2*pi)/3), 0, 0, 0, 0, 0, -lst*cos(phi3)*cos(shi3), 0, 0, 0, 0, 0, lst*sin(phi3)*sin(shi3), 0, 0, 0, 0, 0, 0, 0, - rt*cos(c)*sin(b)*cos(gmt - (2*pi)/3) - rt*sin(b)*sin(c)*sin(gmt - (2*pi)/3), rt*cos(b)*cos(c)*sin(gmt - (2*pi)/3) - rt*cos(b)*sin(c)*cos(gmt - (2*pi)/3); 0, 0, lc*sin(th3)*cos(gmb - (2*pi)/3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, lst*cos(shi3), 0, 0, 0, 0, 0, 0, rt*sin(gmt - (2*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - rt*cos(gmt - (2*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)), rt*cos(b)*cos(c)*sin(a)*cos(gmt - (2*pi)/3) + rt*cos(b)*sin(a)*sin(c)*sin(gmt - (2*pi)/3), rt*cos(gmt - (2*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) + rt*sin(gmt - (2*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)); 0, 0, -lc*cos(th3), 0, 0, 0, 0, 0, lst*cos(shi3)*sin(phi3), 0, 0, 0, 0, 0, lst*cos(phi3)*sin(shi3), 0, 0, 0, 0, 0, 0, rt*cos(gmt - (2*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) - rt*sin(gmt - (2*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)), - rt*cos(a)*cos(b)*cos(c)*cos(gmt - (2*pi)/3) - rt*cos(a)*cos(b)*sin(c)*sin(gmt - (2*pi)/3), rt*cos(gmt - (2*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) + rt*sin(gmt - (2*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)); 0, 0, 0, -lc*sin(th4)*sin(gmb + (2*pi)/3), 0, 0, 0, 0, 0, -lst*cos(phi4)*cos(shi4), 0, 0, 0, 0, 0, lst*sin(phi4)*sin(shi4), 0, 0, 0, 0, 0, 0, rt*sin(b)*sin(c)*sin(gmt + (2*pi)/3) - rt*cos(c)*sin(b)*cos(gmt + (2*pi)/3), - rt*cos(b)*cos(c)*sin(gmt + (2*pi)/3) - rt*cos(b)*sin(c)*cos(gmt + (2*pi)/3); 0, 0, 0, lc*sin(th4)*cos(gmb + (2*pi)/3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, lst*cos(shi4), 0, 0, 0, 0, 0, - rt*cos(gmt + (2*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) - rt*sin(gmt + (2*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)), rt*cos(b)*cos(c)*sin(a)*cos(gmt + (2*pi)/3) - rt*cos(b)*sin(a)*sin(c)*sin(gmt + (2*pi)/3), rt*cos(gmt + (2*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - rt*sin(gmt + (2*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)); 0, 0, 0, -lc*cos(th4), 0, 0, 0, 0, 0, lst*cos(shi4)*sin(phi4), 0, 0, 0, 0, 0, lst*cos(phi4)*sin(shi4), 0, 0, 0, 0, 0, rt*cos(gmt + (2*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) + rt*sin(gmt + (2*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)), rt*cos(a)*cos(b)*sin(c)*sin(gmt + (2*pi)/3) - rt*cos(a)*cos(b)*cos(c)*cos(gmt + (2*pi)/3), rt*cos(gmt + (2*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - rt*sin(gmt + (2*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)); 0, 0, 0, 0, lc*sin(th5)*sin(gmb - (4*pi)/3), 0, 0, 0, 0, 0, -lst*cos(phi5)*cos(shi5), 0, 0, 0, 0, 0, lst*sin(phi5)*sin(shi5), 0, 0, 0, 0, 0, - rt*cos(c)*sin(b)*cos(gmt - (4*pi)/3) - rt*sin(b)*sin(c)*sin(gmt - (4*pi)/3), rt*cos(b)*cos(c)*sin(gmt - (4*pi)/3) - rt*cos(b)*sin(c)*cos(gmt - (4*pi)/3); 0, 0, 0, 0, lc*sin(th5)*cos(gmb - (4*pi)/3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, lst*cos(shi5), 0, 0, 0, 0, rt*sin(gmt - (4*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - rt*cos(gmt - (4*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)), rt*cos(b)*cos(c)*sin(a)*cos(gmt - (4*pi)/3) + rt*cos(b)*sin(a)*sin(c)*sin(gmt - (4*pi)/3), rt*cos(gmt - (4*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) + rt*sin(gmt - (4*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)); 0, 0, 0, 0, -lc*cos(th5), 0, 0, 0, 0, 0, lst*cos(shi5)*sin(phi5), 0, 0, 0, 0, 0, lst*cos(phi5)*sin(shi5), 0, 0, 0, 0, rt*cos(gmt - (4*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) - rt*sin(gmt - (4*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)), - rt*cos(a)*cos(b)*cos(c)*cos(gmt - (4*pi)/3) - rt*cos(a)*cos(b)*sin(c)*sin(gmt - (4*pi)/3), rt*cos(gmt - (4*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) + rt*sin(gmt - (4*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)); 0, 0, 0, 0, 0, -lc*sin(th6)*sin(gmb + (4*pi)/3), 0, 0, 0, 0, 0, -lst*cos(phi6)*cos(shi6), 0, 0, 0, 0, 0, lst*sin(phi6)*sin(shi6), 0, 0, 0, 0, rt*sin(b)*sin(c)*sin(gmt + (4*pi)/3) - rt*cos(c)*sin(b)*cos(gmt + (4*pi)/3), - rt*cos(b)*cos(c)*sin(gmt + (4*pi)/3) - rt*cos(b)*sin(c)*cos(gmt + (4*pi)/3); 0, 0, 0, 0, 0, lc*sin(th6)*cos(gmb + (4*pi)/3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, lst*cos(shi6), 0, 0, 0, - rt*cos(gmt + (4*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) - rt*sin(gmt + (4*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)), rt*cos(b)*cos(c)*sin(a)*cos(gmt + (4*pi)/3) - rt*cos(b)*sin(a)*sin(c)*sin(gmt + (4*pi)/3), rt*cos(gmt + (4*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - rt*sin(gmt + (4*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)); 0, 0, 0, 0, 0, -lc*cos(th6), 0, 0, 0, 0, 0, lst*cos(shi6)*sin(phi6), 0, 0, 0, 0, 0, lst*cos(phi6)*sin(shi6), 0, 0, 0, rt*cos(gmt + (4*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) + rt*sin(gmt + (4*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)), rt*cos(a)*cos(b)*sin(c)*sin(gmt + (4*pi)/3) - rt*cos(a)*cos(b)*cos(c)*cos(gmt + (4*pi)/3), rt*cos(gmt + (4*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - rt*sin(gmt + (4*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b))]


%Aprime=jac*inv(mass)*transpose(jac)
