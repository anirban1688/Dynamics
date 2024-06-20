function TE1=rss_TE(geom,sol,dim,inertia)
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

TE1=zeros(size(sol, 1),1);
varE=zeros(size(sol, 1),1);

for i = 1:size(sol, 1) % Iterate over rows
    % Extract the values of q and q dot from the current row
    q = sol(i, 1:24)
    dq = sol(i, 25:48)

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
        % 
        vth1= dq(1);
        vth2= dq(2);
        vth3= dq(3);
        vth4= dq(4);
        vth5= dq(5);
        vth6= dq(6);

        vphi1= dq(7);
        vphi2= dq(8);
        vphi3= dq(9);
        vphi4= dq(10);
        vphi5= dq(11);
        vphi6= dq(12);

        vshi1= dq(13);
        vshi2= dq(14);
        vshi3= dq(15);
        vshi4= dq(16);
        vshi5= dq(17);
        vshi6= dq(18);

        vx= dq(19);
        vy= dq(20);
        vz= dq(21);
        va= dq(22);
        vb= dq(23);
        vc= dq(24);

        %TE1(i)=49*mst*lc*sin(th1)+ lst*vx*vy*vz*cos(phi1)*cos(shi1);

TE1(i)=(49*mst*(lc*sin(th1) + (lst*vx*vy*vz*cos(phi1)*cos(shi1))/2))/5 + (49*mst*(lc*sin(th2) + (lst*vx*vy*vz*cos(phi2)*cos(shi2))/2))/5 + (49*mst*(lc*sin(th3) + (lst*vx*vy*vz*cos(phi3)*cos(shi3))/2))/5 + (49*mst*(lc*sin(th4) + (lst*vx*vy*vz*cos(phi4)*cos(shi4))/2))/5 + (49*mst*(lc*sin(th5) + (lst*vx*vy*vz*cos(phi5)*cos(shi5))/2))/5 + (49*mst*(lc*sin(th6) + (lst*vx*vy*vz*cos(phi6)*cos(shi6))/2))/5 + (49*mp*z)/5 + Iyy*((vx*vy*vz*(va*cos(b)*sin(a) + vb*cos(a)*sin(b))*(sin(a)*sin(c) - vx*vy*vz*cos(a)*cos(c)*sin(b)))/2 + (vx*vy*vz*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b))*(va*cos(a)*cos(b) - vb*vx*vy*vz*sin(a)*sin(b)))/2 - (vb*vx*vy*vz*cos(b)^2*cos(c))/2)*((va*cos(b)*sin(a) + vb*cos(a)*sin(b))*(sin(a)*sin(c) - vx*vy*vz*cos(a)*cos(c)*sin(b)) + (cos(a)*sin(c) + cos(c)*sin(a)*sin(b))*(va*cos(a)*cos(b) - vb*vx*vy*vz*sin(a)*sin(b)) - vb*vx*vy*vz*cos(b)^2*cos(c)) + Izz*((cos(c)*sin(a) + cos(a)*sin(b)*sin(c))*(va*cos(a)*sin(c) + vc*cos(c)*sin(a) + va*cos(c)*sin(a)*sin(b) + vc*cos(a)*sin(b)*sin(c) - vb*vx*vy*vz*cos(a)*cos(b)*cos(c)) + (cos(a)*cos(c) - vx*vy*vz*sin(a)*sin(b)*sin(c))*(vc*cos(a)*cos(c) + va*cos(a)*cos(c)*sin(b) + vb*cos(b)*cos(c)*sin(a) - va*vx*vy*vz*sin(a)*sin(c) - vc*vx*vy*vz*sin(a)*sin(b)*sin(c)) + cos(b)*sin(c)*(vb*cos(c)*sin(b) + vc*cos(b)*sin(c)))*((vx*vy*vz*(cos(a)*cos(c) - vx*vy*vz*sin(a)*sin(b)*sin(c))*(vc*cos(a)*cos(c) + va*cos(a)*cos(c)*sin(b) + vb*cos(b)*cos(c)*sin(a) - va*vx*vy*vz*sin(a)*sin(c) - vc*vx*vy*vz*sin(a)*sin(b)*sin(c)))/2 + (vx*vy*vz*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c))*(va*cos(a)*sin(c) + vc*cos(c)*sin(a) + va*cos(c)*sin(a)*sin(b) + vc*cos(a)*sin(b)*sin(c) - vb*vx*vy*vz*cos(a)*cos(b)*cos(c)))/2 + (vx*vy*vz*cos(b)*sin(c)*(vb*cos(c)*sin(b) + vc*cos(b)*sin(c)))/2) + Ixx2*vshi1*((vx*vy*vz*(sin(shi1)*cos(gmb) + cos(shi1)*sin(phi1)*sin(gmb))*(vshi1*sin(shi1)*cos(gmb) + vphi1*cos(phi1)*sin(shi1)*sin(gmb) + vshi1*cos(shi1)*sin(phi1)*sin(gmb)))/2 - (vx*vy*vz*(sin(shi1)*sin(gmb) - vx*vy*vz*cos(shi1)*sin(phi1)*cos(gmb))*(vphi1*cos(phi1)*sin(shi1)*cos(gmb) + vshi1*cos(shi1)*sin(phi1)*cos(gmb) - vshi1*vx*vy*vz*sin(shi1)*sin(gmb)))/2 + (vx*vy*vz*cos(phi1)*cos(shi1)*(vshi1*cos(phi1)*cos(shi1) - vphi1*vx*vy*vz*sin(phi1)*sin(shi1)))/2) + Ixx2*vshi2*((vx*vy*vz*(sin(shi2)*sin(gmb) + cos(shi2)*sin(phi2)*cos(gmb))*(vshi2*sin(shi2)*sin(gmb) + vphi2*cos(phi2)*sin(shi2)*cos(gmb) + vshi2*cos(shi2)*sin(phi2)*cos(gmb)))/2 - (vx*vy*vz*(sin(shi2)*cos(gmb) - vx*vy*vz*cos(shi2)*sin(phi2)*sin(gmb))*(vphi2*cos(phi2)*sin(shi2)*sin(gmb) + vshi2*cos(shi2)*sin(phi2)*sin(gmb) - vshi2*vx*vy*vz*sin(shi2)*cos(gmb)))/2 + (vx*vy*vz*cos(phi2)*cos(shi2)*(vshi2*cos(phi2)*cos(shi2) - vphi2*vx*vy*vz*sin(phi2)*sin(shi2)))/2) + (49*lc*mck*sin(th1))/10 + (49*lc*mck*sin(th2))/10 + (49*lc*mck*sin(th3))/10 + (49*lc*mck*sin(th4))/10 + (49*lc*mck*sin(th5))/10 + (49*lc*mck*sin(th6))/10 + Ixx*((vx*vy*vz*cos(b)*sin(a)*(va*cos(c)*sin(a) + vc*cos(a)*sin(c) + va*cos(a)*sin(b)*sin(c) + vb*cos(b)*sin(a)*sin(c) + vc*cos(c)*sin(a)*sin(b)))/2 - (vx*vy*vz*sin(b)*(vc*cos(b)*cos(c) - vb*vx*vy*vz*sin(b)*sin(c)))/2 + (vx*vy*vz*cos(a)*cos(b)*(va*cos(a)*cos(c) + vb*cos(a)*cos(b)*sin(c) + vc*cos(a)*cos(c)*sin(b) - vc*vx*vy*vz*sin(a)*sin(c) - va*vx*vy*vz*sin(a)*sin(b)*sin(c)))/2)*(cos(a)*cos(b)*(va*cos(a)*cos(c) + vb*cos(a)*cos(b)*sin(c) + vc*cos(a)*cos(c)*sin(b) - vc*vx*vy*vz*sin(a)*sin(c) - va*vx*vy*vz*sin(a)*sin(b)*sin(c)) + cos(b)*sin(a)*(va*cos(c)*sin(a) + vc*cos(a)*sin(c) + va*cos(a)*sin(b)*sin(c) + vb*cos(b)*sin(a)*sin(c) + vc*cos(c)*sin(a)*sin(b)) - vx*vy*vz*sin(b)*(vc*cos(b)*cos(c) - vb*vx*vy*vz*sin(b)*sin(c))) + (mst*vx*vy*vz*((vx*vy*vz*(2*lc*vth1*sin(th1)*sin(gmb) + lst*vshi1*sin(phi1)*sin(shi1) - lst*vphi1*vx*vy*vz*cos(phi1)*cos(shi1))^2)/4 + (vx*vy*vz*(lst*vshi1*cos(shi1) + 2*lc*vth1*sin(th1)*cos(gmb))^2)/4 + (vx*vy*vz*(lst*vphi1*cos(shi1)*sin(phi1) - 2*lc*vth1*cos(th1) + lst*vshi1*cos(phi1)*sin(shi1))^2)/4))/2 + (mst*vx*vy*vz*((vx*vy*vz*(2*lc*vth2*sin(th2)*sin(gmb) + lst*vphi2*cos(phi2)*cos(shi2) - lst*vshi2*vx*vy*vz*sin(phi2)*sin(shi2))^2)/4 + (vx*vy*vz*(lst*vshi2*cos(shi2) + 2*lc*vth2*sin(th2)*cos(gmb))^2)/4 + (vx*vy*vz*(lst*vphi2*cos(shi2)*sin(phi2) - 2*lc*vth2*cos(th2) + lst*vshi2*cos(phi2)*sin(shi2))^2)/4))/2 + (vth1^2*vx*vy*vz*(4*Izz1 + lc^2*mck))/8 + (vth2^2*vx*vy*vz*(4*Izz1 + lc^2*mck))/8 + (vth3^2*vx*vy*vz*(4*Izz1 + lc^2*mck))/8 + (vth4^2*vx*vy*vz*(4*Izz1 + lc^2*mck))/8 + (mp*vx*vy*vz*vx^2)/2 + (mp*vx*vy*vz*vy^2)/2 + (mp*vx*vy*vz*vz^2)/2 + (mck*vx*vy*vz*((lc^2*vth5^2*vx*vy*vz*cos(th5)^2)/4 + (lc^2*vth5^2*vx*vy*vz*sin(th5)^2*cos(gmb - (4*pi)/3)^2)/4 + (lc^2*vth5^2*vx*vy*vz*sin(th5)^2*sin(gmb - (4*pi)/3)^2)/4))/2 + (mck*vx*vy*vz*((lc^2*vth6^2*vx*vy*vz*cos(th6)^2)/4 + (lc^2*vth6^2*vx*vy*vz*sin(th6)^2*cos(gmb + (4*pi)/3)^2)/4 + (lc^2*vth6^2*vx*vy*vz*sin(th6)^2*sin(gmb + (4*pi)/3)^2)/4))/2 + (Izz2*vphi3^2*vx*vy*vz)/2 + (Izz2*vphi4^2*vx*vy*vz)/2 + (Izz2*vphi5^2*vx*vy*vz)/2 + (Izz2*vphi6^2*vx*vy*vz)/2 + (Ixx2*vshi3^2*vx*vy*vz)/2 + (Ixx2*vshi4^2*vx*vy*vz)/2 + (Ixx2*vshi5^2*vx*vy*vz)/2 + (Ixx2*vshi6^2*vx*vy*vz)/2 + (Izz1*vth5^2*vx*vy*vz)/2 + (Izz1*vth6^2*vx*vy*vz)/2 + (Iyy2*vphi1^2*vx*vy*vz*cos(shi1)^2)/2 + (Iyy2*vphi2^2*vx*vy*vz*cos(shi2)^2)/2 + (Iyy2*vphi3^2*vx*vy*vz*cos(shi3)^2)/2 + (Iyy2*vphi4^2*vx*vy*vz*cos(shi4)^2)/2 + (Iyy2*vphi5^2*vx*vy*vz*cos(shi5)^2)/2 + (Iyy2*vphi6^2*vx*vy*vz*cos(shi6)^2)/2 - (Izz2*vphi3^2*vx*vy*vz*cos(shi3)^2)/2 - (Izz2*vphi4^2*vx*vy*vz*cos(shi4)^2)/2 - (Izz2*vphi5^2*vx*vy*vz*cos(shi5)^2)/2 - (Izz2*vphi6^2*vx*vy*vz*cos(shi6)^2)/2 + (Izz2*vphi1^2*vx*vy*vz*sin(shi1)^2)/2 + (Izz2*vphi2^2*vx*vy*vz*sin(shi2)^2)/2 + (lc^2*mst*vth3^2*vx*vy*vz)/2 + (lc^2*mst*vth4^2*vx*vy*vz)/2 + (lc^2*mst*vth5^2*vx*vy*vz)/2 + (lc^2*mst*vth6^2*vx*vy*vz)/2 + (lst^2*mst*vshi3^2*vx*vy*vz)/8 + (lst^2*mst*vshi4^2*vx*vy*vz)/8 + (lst^2*mst*vshi5^2*vx*vy*vz)/8 + (lst^2*mst*vshi6^2*vx*vy*vz)/8 + (lst^2*mst*vphi3^2*vx*vy*vz*cos(shi3)^2)/8 + (lst^2*mst*vphi4^2*vx*vy*vz*cos(shi4)^2)/8 + (lst^2*mst*vphi5^2*vx*vy*vz*cos(shi5)^2)/8 + (lst^2*mst*vphi6^2*vx*vy*vz*cos(shi6)^2)/8 - (lc*lst*mst*vshi3*vth3*vx*vy*vz*cos(shi3)*sin(th3)*cos(gmb))/4 - (lc*lst*mst*vshi4*vth4*vx*vy*vz*cos(shi4)*sin(th4)*cos(gmb))/4 - (lc*lst*mst*vshi5*vth5*vx*vy*vz*cos(shi5)*sin(th5)*cos(gmb))/4 - (lc*lst*mst*vshi6*vth6*vx*vy*vz*cos(shi6)*sin(th6)*cos(gmb))/4 - (lc*lst*mst*vphi3*vth3*vx*vy*vz*cos(shi3)*cos(th3)*sin(phi3))/2 - (lc*lst*mst*vphi4*vth4*vx*vy*vz*cos(shi4)*cos(th4)*sin(phi4))/2 - (lc*lst*mst*vphi5*vth5*vx*vy*vz*cos(shi5)*cos(th5)*sin(phi5))/2 - (lc*lst*mst*vphi6*vth6*vx*vy*vz*cos(shi6)*cos(th6)*sin(phi6))/2 - (lc*lst*mst*vshi3*vth3*vx*vy*vz*cos(phi3)*cos(th3)*sin(shi3))/2 - (lc*lst*mst*vshi4*vth4*vx*vy*vz*cos(phi4)*cos(th4)*sin(shi4))/2 - (lc*lst*mst*vshi5*vth5*vx*vy*vz*cos(phi5)*cos(th5)*sin(shi5))/2 - (lc*lst*mst*vshi6*vth6*vx*vy*vz*cos(phi6)*cos(th6)*sin(shi6))/2 + (3^((vx*vy*vz)/2)*lc*lst*mst*vshi3*vth3*vx*vy*vz*cos(shi3)*sin(th3)*sin(gmb))/4 - (3^((vx*vy*vz)/2)*lc*lst*mst*vshi4*vth4*vx*vy*vz*cos(shi4)*sin(th4)*sin(gmb))/4 - (3^((vx*vy*vz)/2)*lc*lst*mst*vshi5*vth5*vx*vy*vz*cos(shi5)*sin(th5)*sin(gmb))/4 + (3^((vx*vy*vz)/2)*lc*lst*mst*vshi6*vth6*vx*vy*vz*cos(shi6)*sin(th6)*sin(gmb))/4 + (lc*lst*mst*vphi3*vth3*vx*vy*vz*cos(phi3)*cos(shi3)*sin(th3)*sin(gmb))/4 - (lc*lst*mst*vphi4*vth4*vx*vy*vz*cos(phi4)*cos(shi4)*sin(th4)*sin(gmb))/4 + (lc*lst*mst*vphi5*vth5*vx*vy*vz*cos(phi5)*cos(shi5)*sin(th5)*sin(gmb))/4 - (lc*lst*mst*vphi6*vth6*vx*vy*vz*cos(phi6)*cos(shi6)*sin(th6)*sin(gmb))/4 - (lc*lst*mst*vshi3*vth3*vx*vy*vz*sin(phi3)*sin(shi3)*sin(th3)*sin(gmb))/4 + (lc*lst*mst*vshi4*vth4*vx*vy*vz*sin(phi4)*sin(shi4)*sin(th4)*sin(gmb))/4 - (lc*lst*mst*vshi5*vth5*vx*vy*vz*sin(phi5)*sin(shi5)*sin(th5)*sin(gmb))/4 + (lc*lst*mst*vshi6*vth6*vx*vy*vz*sin(phi6)*sin(shi6)*sin(th6)*sin(gmb))/4 + (3^((vx*vy*vz)/2)*lc*lst*mst*vphi3*vth3*vx*vy*vz*cos(phi3)*cos(shi3)*sin(th3)*cos(gmb))/4 + (3^((vx*vy*vz)/2)*lc*lst*mst*vphi4*vth4*vx*vy*vz*cos(phi4)*cos(shi4)*sin(th4)*cos(gmb))/4 - (3^((vx*vy*vz)/2)*lc*lst*mst*vphi5*vth5*vx*vy*vz*cos(phi5)*cos(shi5)*sin(th5)*cos(gmb))/4 - (3^((vx*vy*vz)/2)*lc*lst*mst*vphi6*vth6*vx*vy*vz*cos(phi6)*cos(shi6)*sin(th6)*cos(gmb))/4 - (3^((vx*vy*vz)/2)*lc*lst*mst*vshi3*vth3*vx*vy*vz*sin(phi3)*sin(shi3)*sin(th3)*cos(gmb))/4 - (3^((vx*vy*vz)/2)*lc*lst*mst*vshi4*vth4*vx*vy*vz*sin(phi4)*sin(shi4)*sin(th4)*cos(gmb))/4 + (3^((vx*vy*vz)/2)*lc*lst*mst*vshi5*vth5*vx*vy*vz*sin(phi5)*sin(shi5)*sin(th5)*cos(gmb))/4 + (3^((vx*vy*vz)/2)*lc*lst*mst*vshi6*vth6*vx*vy*vz*sin(phi6)*sin(shi6)*sin(th6)*cos(gmb))/4;
end
%plot(TE1)
% for i = 1:size(sol, 1)
%     varE(i)=(TE1(1)-TE1(i))/TE1(1);
% end
% plot(TE1)