function dydt=diffeq3(t,y)
dydt=zeros(48,1);
 
    % Extract the values of q and q dot from the current row
    q =y(1:24)
    dq =y(25:48)

rb = 0.22774;   % Fixed platform radius in m
rt = 0.186;   % Moving platform radius in m
gmb = 0.406; % Base angle in radians
gmt = 0.109; % Top angle in radians
lst = 0.182;  % Length of the struts in m
lc = 0.0125;    % Length of the cranks in m
geom=[rb;rt;gmb;gmt;lst;lc];


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

    M=rss_mass2(geom,inertia,q,dim);
    C=rss_cor2(geom,inertia,q,dq,dim);
    J=rss_jac2(geom,q,dim);
    A=J*inv(M)*transpose(J);
    G=rss_G1(geom,inertia,q,dim);
    Jdot=rss_jacdot2(geom,q,dq,dim);


    Ainv=inv(A);
    F=-C*dq-transpose(G);
    Lmd=-Ainv*(Jdot*dq+J*inv(M)*F);
    qddot=inv(M)*F+inv(M)*transpose(J)*Lmd;



    vector1=[y(25) ; y(26) ; y(27) ; y(28) ; y(29) ; y(30) ; y(31) ; y(32) ; y(33) ; y(34) ; y(35) ; y(36) ; y(37) ; y(38) ; y(39) ; y(40) ; y(41) ; y(42) ; y(43) ; y(44) ; y(45) ; y(46) ; y(47) ; y(48)]
    vector2=[qddot];
    dydt=[vector1;vector2];
end
