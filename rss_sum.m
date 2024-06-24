
function rss_sum(Q,geom)



rb = geom(1);
rt = geom(2);
gmb= geom(3);
gmt = geom(4);
lst = geom(5);
lc = geom(6);
sum=zeros(size(Q, 1),1);


for i = 1:size(Q, 1) % Iterate over rows
    % Extract the values of q and q dot from the current row
    q = Q(i, 1:24)
    dq = Q(i, 25:48)

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

      
        ceq=[x - rb*cos(gmb) - lst*cos(shi1)*sin(phi1) - lc*cos(th1)*sin(gmb) + rt*cos(b)*cos(c)*cos(gmt) + rt*cos(b)*sin(c)*sin(gmt); y + lst*sin(shi1) + rb*sin(gmb) + rt*cos(gmt)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) - rt*sin(gmt)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - lc*cos(th1)*cos(gmb); z - lc*sin(th1) + rt*cos(gmt)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) - rt*sin(gmt)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - lst*cos(phi1)*cos(shi1); x - rb*cos(gmb) - lst*cos(shi2)*sin(phi2) + lc*cos(th2)*sin(gmb) + rt*cos(b)*cos(c)*cos(gmt) - rt*cos(b)*sin(c)*sin(gmt); y + lst*sin(shi2) - rb*sin(gmb) + rt*cos(gmt)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) + rt*sin(gmt)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - lc*cos(th2)*cos(gmb); z - lc*sin(th2) + rt*cos(gmt)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) + rt*sin(gmt)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - lst*cos(phi2)*cos(shi2); x - rb*cos(gmb - (2*pi)/3) - lc*cos(th3)*sin(gmb - (2*pi)/3) - lst*cos(shi3)*sin(phi3) + rt*cos(b)*cos(c)*cos(gmt - (2*pi)/3) + rt*cos(b)*sin(c)*sin(gmt - (2*pi)/3); y + lst*sin(shi3) + rb*sin(gmb - (2*pi)/3) + rt*cos(gmt - (2*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) - rt*sin(gmt - (2*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - lc*cos(th3)*cos(gmb - (2*pi)/3); z - lc*sin(th3) + rt*cos(gmt - (2*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) - rt*sin(gmt - (2*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - lst*cos(phi3)*cos(shi3); x - rb*cos(gmb + (2*pi)/3) + lc*cos(th4)*sin(gmb + (2*pi)/3) - lst*cos(shi4)*sin(phi4) + rt*cos(b)*cos(c)*cos(gmt + (2*pi)/3) - rt*cos(b)*sin(c)*sin(gmt + (2*pi)/3); y + lst*sin(shi4) - rb*sin(gmb + (2*pi)/3) + rt*cos(gmt + (2*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) + rt*sin(gmt + (2*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - lc*cos(th4)*cos(gmb + (2*pi)/3); z - lc*sin(th4) + rt*cos(gmt + (2*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) + rt*sin(gmt + (2*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - lst*cos(phi4)*cos(shi4); x - rb*cos(gmb - (4*pi)/3) - lc*cos(th5)*sin(gmb - (4*pi)/3) - lst*cos(shi5)*sin(phi5) + rt*cos(b)*cos(c)*cos(gmt - (4*pi)/3) + rt*cos(b)*sin(c)*sin(gmt - (4*pi)/3); y + lst*sin(shi5) + rb*sin(gmb - (4*pi)/3) + rt*cos(gmt - (4*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) - rt*sin(gmt - (4*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - lc*cos(th5)*cos(gmb - (4*pi)/3); z - lc*sin(th5) + rt*cos(gmt - (4*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) - rt*sin(gmt - (4*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - lst*cos(phi5)*cos(shi5); x - rb*cos(gmb + (4*pi)/3) + lc*cos(th6)*sin(gmb + (4*pi)/3) - lst*cos(shi6)*sin(phi6) + rt*cos(b)*cos(c)*cos(gmt + (4*pi)/3) - rt*cos(b)*sin(c)*sin(gmt + (4*pi)/3); y + lst*sin(shi6) - rb*sin(gmb + (4*pi)/3) + rt*cos(gmt + (4*pi)/3)*(cos(a)*sin(c) + cos(c)*sin(a)*sin(b)) + rt*sin(gmt + (4*pi)/3)*(cos(a)*cos(c) - sin(a)*sin(b)*sin(c)) - lc*cos(th6)*cos(gmb + (4*pi)/3); z - lc*sin(th6) + rt*cos(gmt + (4*pi)/3)*(sin(a)*sin(c) - cos(a)*cos(c)*sin(b)) + rt*sin(gmt + (4*pi)/3)*(cos(c)*sin(a) + cos(a)*sin(b)*sin(c)) - lst*cos(phi6)*cos(shi6)]
        sum(i)=transpose(ceq)*ceq;


end

plot(sum)
