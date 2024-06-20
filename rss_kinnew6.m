 clc
clear all

% % % Platform and Joint Parameters
rb_value = 227.74;   % Fixed platform radius in mm
rt_value = 186;   % Moving platform radius in mm
gmb_value = 0.406; % Base angle in radians
gmt_value = 0.109; % Top angle in radians
z_value = 160;    % Height of the moving platform from base in mm
x_value = 0;    % X-position of the platform in mm
y_value = 0;     % Y-position of the platform in mm
lst_value = 182;  % Length of the struts in mm
lc_value = 12.5;    % Length of the cranks in mm

% % Rotation Parameters
a_value = 0;      % Rotation about x-axis in radians
b_value = 0;      % Rotation about y-axis in radians
c_value = 0;      % Rotation about z-axis in radians

syms gmt gmb x y z a b c rt rb th1 th2 th3 th4 th5 th6 lc lst ta1 ta2 ta3 ta4 ta5 ta6    
syms phi1 phi2 phi3 phi4 phi5 phi6 shi1 shi2 shi3 shi4 shi5 shi6
Rxa=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)]
Ryb=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)]
Rzc=[cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1]
Rxyz= Rxa*Ryb*Rzc;

gmb1=-gmb;
gmb2=gmb;
gmb3=2*pi/3-gmb;
gmb4=2*pi/3+gmb;
gmb5=4*pi/3-gmb;
gmb6=4*pi/3+gmb;


%%vertices of the fixed platform(bi)w.r.t global frame of reference%%
b1=rb*[cos(-gmb);sin(-gmb);0]
b2=rb*[cos(gmb);sin(gmb);0]
b3=rb*[cos(2*pi/3-gmb);sin(2*pi/3-gmb);0]
b4=rb*[cos(2*pi/3+gmb);sin(2*pi/3+gmb);0]
b5=rb*[cos(4*pi/3-gmb);sin(4*pi/3-gmb);0]
b6=rb*[cos(4*pi/3+gmb);sin(4*pi/3+gmb);0]



%%vertices of the moving platform(ti)w.r.t local frame of reference%%
t1=rt*[cos(-gmt);sin(-gmt);0]
t2=rt*[cos(gmt);sin(gmt);0]
t3=rt*[cos(2*pi/3-gmt);sin(2*pi/3-gmt);0]
t4=rt*[cos(2*pi/3+gmt);sin(2*pi/3+gmt);0]
t5=rt*[cos(4*pi/3-gmt);sin(4*pi/3-gmt);0]
t6=rt*[cos(4*pi/3+gmt);sin(4*pi/3+gmt);0]

t=[t1 t2 t3 t4 t5 t6];

%%tiprime is the new co-ordinate of vertices of the moving platform
%% w.r.t global frame of reference
P=[x;y;z]
for j= 1:6
tiprime=P+Rxyz*t(:,j)
end

%gmbvec=[gmb1 gmb2 gmb3 gmb4 gmb5 gmb6]
%%rotation matrix about z %%
    Rz1=[cos(gmb1) -sin(gmb1) 0 ; sin(gmb1) cos(gmb1) 0 ; 0 0 1]
    Rz2=[cos(gmb2) -sin(gmb2) 0;sin(gmb2) cos(gmb2) 0;0 0 1]
    Rz3=[cos(gmb3) -sin(gmb3) 0;sin(gmb3) cos(gmb3) 0;0 0 1]
    Rz4=[cos(gmb4) -sin(gmb4) 0;sin(gmb4) cos(gmb4) 0;0 0 1]
    Rz5=[cos(gmb5) -sin(gmb5) 0;sin(gmb5) cos(gmb5) 0;0 0 1]
    Rz6=[cos(gmb6) -sin(gmb6) 0;sin(gmb6) cos(gmb6) 0;0 0 1]




Rz=[Rz1 Rz2 Rz3 Rz4 Rz5 Rz6];
th=[th1 th2 th3 th4 th5 th6];

%%centre of mass of Cranks%%
ck1=Rz1*[rb;(lc)*cos(th1);(lc)*sin(th1)]
ck2=Rz2*[rb;(lc)*cos(th2);(lc)*sin(th2)]
ck3=Rz3*[rb;(lc)*cos(th3);(lc)*sin(th3)]
ck4=Rz4*[rb;(lc)*cos(th4);(lc)*sin(th4)]
ck5=Rz5*[rb;(lc)*cos(th5);(lc)*sin(th5)]
ck6=Rz6*[rb;(lc)*cos(th6);(lc)*sin(th6)]

%%rotation matrix aboutx %%
Rx1=[1 0 0;0 cos(shi1) -sin(shi1);0 sin(shi1) cos(shi1)]
Rx2=[1 0 0;0 cos(shi2) -sin(shi2);0 sin(shi2) cos(shi2)]
Rx3=[1 0 0;0 cos(shi3) -sin(shi3);0 sin(shi3) cos(shi3)]
Rx4=[1 0 0;0 cos(shi4) -sin(shi4);0 sin(shi4) cos(shi4)]
Rx5=[1 0 0;0 cos(shi5) -sin(shi5);0 sin(shi5) cos(shi5)]
Rx6=[1 0 0;0 cos(shi6) -sin(shi6);0 sin(shi6) cos(shi6)]

%%rotation matrix about y %%
Ry1=[cos(phi1) 0 sin(phi1);0 1 0;-sin(phi1) 0 cos(phi1)]
Ry2=[cos(phi2) 0 sin(phi2);0 1 0;-sin(phi2) 0 cos(phi2)]
Ry3=[cos(phi3) 0 sin(phi3);0 1 0;-sin(phi3) 0 cos(phi3)]
Ry4=[cos(phi4) 0 sin(phi4);0 1 0;-sin(phi4) 0 cos(phi4)]
Ry5=[cos(phi5) 0 sin(phi5);0 1 0;-sin(phi5) 0 cos(phi5)]
Ry6=[cos(phi6) 0 sin(phi6);0 1 0;-sin(phi6) 0 cos(phi6)]

Rx=[Rx1 Rx2 Rx3 Rx4 Rx5 Rx6];
Ry=[Ry1 Ry2 Ry3 Ry4 Ry5 Ry6];

cst1=ck1+Ry1*Rx1*[0;0;lst]
cst2=ck2+Ry2*Rx2*[0;0;lst]
cst3=ck3+Ry3*Rx3*[0;0;lst]
cst4=ck4+Ry4*Rx4*[0;0;lst]
cst5=ck5+Ry5*Rx5*[0;0;lst]
cst6=ck6+Ry6*Rx6*[0;0;lst]

cst=[cst1 cst2 cst3 cst4 cst5 cst6]

ceq1=P+Rxyz*t1-cst1
ceq2=P+Rxyz*t2-cst2
ceq3=P+Rxyz*t3-cst3
ceq4=P+Rxyz*t4-cst4
ceq5=P+Rxyz*t5-cst5
ceq6=P+Rxyz*t6-cst6

ceq=[ceq1;ceq2;ceq3;ceq4;ceq5;ceq6]
pos=[th1 th2 th3 th4 th5 th6 phi1 phi2 phi3 phi4 phi5 phi6 shi1 shi2 shi3 shi4 shi5 shi6 x y z a b c]
j=jacobian(ceq,pos)


%% For calculation of theta
vec7=P+Rxyz*t1-ck1;
eq7=simplify(expand(transpose(vec7)*vec7-lst^2));
eq7a=subs(eq7,[cos(th1),sin(th1)],[(1-ta1^2)/(1+ta1^2),2*ta1/(1+ta1^2)]);
solth1=solve(eq7a==0,ta1);
solth1_subs = subs(solth1, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value])
th1_value=double(2*atan(solth1_subs(1)))

vec8=P+Rxyz*t2-ck2;
eq8=simplify(expand(transpose(vec8)*vec8-lst^2));
eq8a=subs(eq8,[cos(th2),sin(th2)],[(1-ta2^2)/(1+ta2^2),2*ta2/(1+ta2^2)]);
solth2=solve(eq8a==0,ta2);
solth2_subs = subs(solth2, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value])
th2_value=double(2*atan(solth2_subs(1)))

vec9=P+Rxyz*t3-ck3;
eq9=simplify(expand(transpose(vec9)*vec9-lst^2));
eq9a=subs(eq9,[cos(th3),sin(th3)],[(1-ta3^2)/(1+ta3^2),2*ta3/(1+ta3^2)]);
solth3=solve(eq9a==0,ta3);
solth3_subs = subs(solth3, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value])
th3_value=double(2*atan(solth3_subs(1)))

vec10=P+Rxyz*t4-ck4;
eq10=simplify(expand(transpose(vec10)*vec10-lst^2));
eq10a=subs(eq10,[cos(th4),sin(th4)],[(1-ta4^2)/(1+ta4^2),2*ta4/(1+ta4^2)]);
solth4=solve(eq10a==0,ta4);
solth4_subs = subs(solth4, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value])
th4_value=double(2*atan(solth4_subs(1)))

vec11=P+Rxyz*t5-ck5;
eq11=simplify(expand(transpose(vec11)*vec11-lst^2));
eq11a=subs(eq11,[cos(th5),sin(th5)],[(1-ta5^2)/(1+ta5^2),2*ta5/(1+ta5^2)]);
solth5=solve(eq11a==0,ta5);
solth5_subs = subs(solth5, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value])
th5_value=double(2*atan(solth5_subs(1)))

vec12=P+Rxyz*t6-ck6;
eq12=simplify(expand(transpose(vec12)*vec12-lst^2));
eq12a=(subs(eq12,[cos(th6),sin(th6)],[(1-ta6^2)/(1+ta6^2),2*ta6/(1+ta6^2)]));
solth6=solve(eq12a==0,ta6);
solth6_subs = subs(solth6, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value])
th6_value=double(2*atan(solth6_subs(1)))

syms s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 
%% for calculation of shi
ceq1=P+Rxyz*t1-cst1;
ceq1a=subs(ceq1(2),sin(shi1),s1);
shi1a=asin(solve(ceq1a==0,s1))
shi1_value=double(subs(shi1a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th1], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th1_value]))

ceq2=P+Rxyz*t2-cst2;
ceq2a=subs(ceq2(2),sin(shi2),s2);
shi2a=asin(solve(ceq2a==0,s2))
shi2_value=double(subs(shi2a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th2], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th2_value]))

ceq3=P+Rxyz*t3-cst3;
ceq3a=subs(ceq3(2),sin(shi3),s3);
shi3a=asin(solve(ceq3a==0,s3))
shi3_value=double(subs(shi3a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th3], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th3_value]))

ceq4=P+Rxyz*t4-cst4;
ceq4a=subs(ceq4(2),sin(shi4),s4);
shi4a=asin(solve(ceq4a==0,s4))
shi4_value=double(subs(shi4a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th4], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th4_value]))

ceq5=P+Rxyz*t5-cst5;
ceq5a=subs(ceq5(2),sin(shi5),s5);
shi5a=asin(solve(ceq5a==0,s5))
shi5_value=double(subs(shi5a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th5], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th5_value]))

ceq6=P+Rxyz*t6-cst6;
ceq6a=subs(ceq6(2),sin(shi6),s6);
shi6a=asin(solve(ceq6a==0,s6))
shi6_value=double(subs(shi6a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th6], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th6_value]))

%% for calculation of phi

ceq1=P+Rxyz*t1-cst1;
ceq7a=subs(ceq1(3),cos(phi1),s7);
phi1a=asin(solve(ceq7a==0,s7))
phi1_value=double(subs(phi1a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th1,shi1], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th1_value,shi1_value]))


ceq2=P+Rxyz*t2-cst2;
ceq8a=subs(ceq2(3),cos(phi2),s8);
phi2a=asin(solve(ceq8a==0,s8))
phi2_value=double(subs(phi2a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th2,shi2], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th2_value,shi2_value]))

ceq3=P+Rxyz*t3-cst3;
ceq9a=subs(ceq3(3),cos(phi3),s9);
phi3a=asin(solve(ceq9a==0,s9))
phi3_value=double(subs(phi3a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th3,shi3], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th3_value,shi3_value]))

ceq4=P+Rxyz*t4-cst4;
ceq10a=subs(ceq4(3),cos(phi4),s10);
phi4a=asin(solve(ceq10a==0,s10))
phi4_value=double(subs(phi4a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th4,shi4], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th4_value,shi4_value]))

ceq5=P+Rxyz*t5-cst5;
ceq11a=subs(ceq5(3),cos(phi5),s11);
phi5a=asin(solve(ceq11a==0,s11))
phi5_value=double(subs(phi5a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th5,shi5], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th5_value,shi5_value]))

ceq6=P+Rxyz*t6-cst6;
ceq12a=subs(ceq6(3),cos(phi6),s12);
phi6a=asin(solve(ceq12a==0,s12))
phi6_value=double(subs(phi6a, [rb, rt, gmb, z, x, y, gmt, lst, lc,a,b,c,th6,shi6], [rb_value, rt_value, gmb_value, z_value, x_value, y_value, gmt_value, lst_value, lc_value,a_value,b_value,c_value,th6_value,shi6_value]))

% q=[th1_value;th2_value;th3_value;th4_value;th5_value;th6_value;shi1_value;shi2_value;shi3_value;shi4_value;shi5_value;shi6_value;phi1_value;phi2_value;phi3_value;phi4_value;phi5_value;phi6_value;x_value;y_value;z_value;a_value;b_value;c    _value]

fid = fopen('solth1.txt', 'wt');
fprintf(fid, '%s\n', char(solth1));
fclose(fid);

fid = fopen('solshi1.txt', 'wt');
fprintf(fid, '%s\n', char(shi1a));
fclose(fid);

fid = fopen('solphi1.txt', 'wt');
fprintf(fid, '%s\n', char(phi1a));
fclose(fid);


fid = fopen('solth2.txt', 'wt');
fprintf(fid, '%s\n', char(solth2));
fclose(fid);

fid = fopen('solshi2.txt', 'wt');
fprintf(fid, '%s\n', char(shi2a));
fclose(fid);

fid = fopen('solphi2.txt', 'wt');
fprintf(fid, '%s\n', char(phi2a));
fclose(fid);

fid = fopen('solth3.txt', 'wt');
fprintf(fid, '%s\n', char(solth3));
fclose(fid);

fid = fopen('solshi3.txt', 'wt');
fprintf(fid, '%s\n', char(shi3a));
fclose(fid);

fid = fopen('solphi3.txt', 'wt');
fprintf(fid, '%s\n', char(phi3a));
fclose(fid);

fid = fopen('solth4.txt', 'wt');
fprintf(fid, '%s\n', char(solth4));
fclose(fid);

fid = fopen('solshi4.txt', 'wt');
fprintf(fid, '%s\n', char(shi4a));
fclose(fid);

fid = fopen('solphi4.txt', 'wt');
fprintf(fid, '%s\n', char(phi4a));
fclose(fid);

fid = fopen('solth5.txt', 'wt');
fprintf(fid, '%s\n', char(solth5));
fclose(fid);

fid = fopen('solshi5.txt', 'wt');
fprintf(fid, '%s\n', char(shi5a));
fclose(fid);

fid = fopen('solphi5.txt', 'wt');
fprintf(fid, '%s\n', char(phi5a));
fclose(fid);


fid = fopen('solth6.txt', 'wt');
fprintf(fid, '%s\n', char(solth6));
fclose(fid);

fid = fopen('solshi6.txt', 'wt');
fprintf(fid, '%s\n', char(shi6a));
fclose(fid);

fid = fopen('solphi6.txt', 'wt');
fprintf(fid, '%s\n', char(phi6a));
fclose(fid);