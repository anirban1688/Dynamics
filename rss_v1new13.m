clc
clear all

syms t a(t) b(t) c(t) th1(t) th2(t) th3(t) th4(t) th5(t) th6(t) phi1(t) phi2(t) phi3(t) phi4(t) phi5(t) phi6(t) shi1(t) shi2(t) shi3(t) shi4(t) shi5(t) shi6(t) 
syms Vap Vbp Vcp lst lc
syms rb rt gmb gmb1 gmb2 gmb3 gmb4 gmb5 gmb6 gmt  x(t) y(t) z(t) A
syms tiprime Rz ckprime cst1dot mck mst
syms t Ixx Iyy Izz Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 omega mp mck mst
syms vth1 vth2 vth3 vth4 vth5 vth6 vphi1 vphi2 vphi3 vphi4 vphi5 vphi6 vshi1 vshi2 vshi3 vshi4 vshi5 vshi6 vx vy vz va vb vc

%% definition of rotation matrix%%%%%%%%%%%
Rxa=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)]
Ryb=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)]
Rzc=[cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1]
Rxyz= Rxa*Ryb*Rzc;

%% angular velocity of the moving platform %%%%%%%%%%%%%
eq=diff(Rxyz,t)
eq1=subs(eq,diff(a(t), t),Vap)
eq2=subs(eq1,diff(b(t), t),Vbp)
dRxyz=subs(eq2,diff(c(t), t),Vcp)
B=simplify(transpose(Rxyz)*dRxyz)    %A=R-dot * transpose of rotation matrix
A=transpose(Rxyz)*dRxyz;
C=simplify(A+transpose(A))         %A is a skew matrix
D=A(t);




%% Angles btw the line (joining each vertices of fixed flatform) and x axis in global frame of reference
gmb1=-gmb;
gmb2=gmb;
gmb3=2*pi/3-gmb;
gmb4=2*pi/3+gmb;
gmb5=4*pi/3-gmb;
gmb6=4*pi/3+gmb;


%% vertices of the fixed platform(bi)w.r.t global frame of reference%%
b1=rb*[cos(-gmb);sin(-gmb);0]
b2=rb*[cos(gmb);sin(gmb);0]
b3=rb*[cos(2*pi/3-gmb);sin(2*pi/3-gmb);0]
b4=rb*[cos(2*pi/3+gmb);sin(2*pi/3+gmb);0]
b5=rb*[cos(4*pi/3-gmb);sin(4*pi/3-gmb);0]
b6=rb*[cos(4*pi/3+gmb);sin(4*pi/3+gmb);0]

%% vertices of the moving platform(ti)w.r.t local frame of reference%%
t1=rt*[cos(-gmt);sin(-gmt);0]
t2=rt*[cos(gmt);sin(gmt);0]
t3=rt*[cos(2*pi/3-gmt);sin(2*pi/3-gmt);0]
t4=rt*[cos(2*pi/3+gmt);sin(2*pi/3+gmt);0]
t5=rt*[cos(4*pi/3-gmt);sin(4*pi/3-gmt);0]
t6=rt*[cos(4*pi/3+gmt);sin(4*pi/3+gmt);0]

t=[t1 t2 t3 t4 t5 t6];

%% tiprime is the new co-ordinate of vertices of the moving platform
%%  w.r.t global frame of reference
P=[x;y;z]
for j= 1:6
tiprime=P+Rxyz*t(:,j)
end





%% rotation matrix about z %%
    Rz1=[cos(gmb1) -sin(gmb1) 0 ; sin(gmb1) cos(gmb1) 0 ; 0 0 1]
    Rz2=[cos(gmb2) -sin(gmb2) 0;sin(gmb2) cos(gmb2) 0;0 0 1]
    Rz3=[cos(gmb3) -sin(gmb3) 0;sin(gmb3) cos(gmb3) 0;0 0 1]
    Rz4=[cos(gmb4) -sin(gmb4) 0;sin(gmb4) cos(gmb4) 0;0 0 1]
    Rz5=[cos(gmb5) -sin(gmb5) 0;sin(gmb5) cos(gmb5) 0;0 0 1]
    Rz6=[cos(gmb6) -sin(gmb6) 0;sin(gmb6) cos(gmb6) 0;0 0 1]




Rz=[Rz1 Rz2 Rz3 Rz4 Rz5 Rz6];
th=[th1 th2 th3 th4 th5 th6];

%% centre of mass of Cranks%%
ck1=Rz1*[rb;(lc/2)*cos(th1);(lc/2)*sin(th1)]
ck2=Rz2*[rb;(lc/2)*cos(th2);(lc/2)*sin(th2)]
ck3=Rz3*[rb;(lc/2)*cos(th3);(lc/2)*sin(th3)]
ck4=Rz4*[rb;(lc/2)*cos(th4);(lc/2)*sin(th4)]
ck5=Rz5*[rb;(lc/2)*cos(th5);(lc/2)*sin(th5)]
ck6=Rz6*[rb;(lc/2)*cos(th6);(lc/2)*sin(th6)]
 
 
  
ck=[ck1 ck2 ck3 ck4 ck5 ck6];

%% velocities of centre of mass of cranks%%
ck1dot=diff(ck1)
ck2dot=diff(ck2)
ck3dot=diff(ck3)
ck4dot=diff(ck4)
ck5dot=diff(ck5)
ck6dot=diff(ck6)

 %% Angular velocity For cranks %%
Wck1=[0;0;diff(th1)]
Wck2=[0;0;diff(th2)]
Wck3=[0;0;diff(th3)]

Wck4=[0;0;diff(th4)]
Wck5=[0;0;diff(th5)]
Wck6=[0;0;diff(th6)]
Wck=[Wck1 Wck2 Wck3 Wck4 Wck5 Wck6]
%% Inertia constant of crank and kinetic energy of each crank
Ic=[Ixx1 0 0;0 Iyy1 0;0 0 Izz1]
KEck1=simplify(0.5*mck*sum(ck1dot.^2) + 0.5*transpose(Wck1)*Ic*Wck1)
simplify(KEck1)
KEck2=simplify(0.5*mck*sum(ck2dot.^2) + 0.5*transpose(Wck2)*Ic*Wck2)
simplify(KEck2)
KEck3=simplify(0.5*mck*sum(ck3dot.^2) + 0.5*transpose(Wck3)*Ic*Wck3)
simplify(KEck3)
KEck4=simplify(0.5*mck*sum(ck4dot.^2) + 0.5*transpose(Wck4)*Ic*Wck4)
simplify(KEck4)
KEck5=0.5*mck*sum(ck5dot.^2) + 0.5*transpose(Wck5)*Ic*Wck5
simplify(KEck5)
KEck6=0.5*mck*sum(ck6dot.^2) + 0.5*transpose(Wck6)*Ic*Wck6
simplify(KEck6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% rotation matrix about x %%
Rx1=[1 0 0;0 cos(shi1) -sin(shi1);0 sin(shi1) cos(shi1)]
Rx2=[1 0 0;0 cos(shi2) -sin(shi2);0 sin(shi2) cos(shi2)]

Rx3=[1 0 0;0 cos(shi3) -sin(shi3);0 sin(shi3) cos(shi3)]

Rx4=[1 0 0;0 cos(shi4) -sin(shi4);0 sin(shi4) cos(shi4)]

Rx5=[1 0 0;0 cos(shi5) -sin(shi5);0 sin(shi5) cos(shi5)]

Rx6=[1 0 0;0 cos(shi6) -sin(shi6);0 sin(shi6) cos(shi6)]

%% rotation matrix about y %%
Ry1=[cos(phi1) 0 sin(phi1);0 1 0;-sin(phi1) 0 cos(phi1)]
Ry2=[cos(phi2) 0 sin(phi2);0 1 0;-sin(phi2) 0 cos(phi2)]
Ry3=[cos(phi3) 0 sin(phi3);0 1 0;-sin(phi3) 0 cos(phi3)]
Ry4=[cos(phi4) 0 sin(phi4);0 1 0;-sin(phi4) 0 cos(phi4)]
Ry5=[cos(phi5) 0 sin(phi5);0 1 0;-sin(phi5) 0 cos(phi5)]
Ry6=[cos(phi6) 0 sin(phi6);0 1 0;-sin(phi6) 0 cos(phi6)]

Rx=[Rx1 Rx2 Rx3 Rx4 Rx5 Rx6];
Ry=[Ry1 Ry2 Ry3 Ry4 Ry5 Ry6];

ck7=Rz1*[rb;lc*cos(th1);lc*sin(th1)]
ck8=Rz2*[rb;(lc)*cos(th2);(lc)*sin(th2)]
ck9=Rz3*[rb;(lc)*cos(th3);(lc)*sin(th3)]
ck10=Rz4*[rb;(lc)*cos(th4);(lc)*sin(th4)]
ck11=Rz5*[rb;(lc)*cos(th5);(lc)*sin(th5)]
ck12=Rz6*[rb;(lc)*cos(th6);(lc)*sin(th6)]
ck13=[ck7 ck8 ck9 ck10 ck11 ck12] 
%% centre of mass of struts
cst1=ck7+Ry1*Rx1*[0;0;lst/2]
cst2=ck8+Ry2*Rx2*[0;0;lst/2]
cst3=ck9+Ry3*Rx3*[0;0;lst/2]
cst4=ck10+Ry4*Rx4*[0;0;lst/2]
cst5=ck11+Ry5*Rx5*[0;0;lst/2]
cst6=ck12+Ry6*Rx6*[0;0;lst/2]

cst=[cst1 cst2 cst3 cst4 cst5 cst6]

%% velocities of centre of mass of struts%%
cst1dot=diff(cst1)
cst2dot=diff(cst2)
cst3dot=diff(cst3)
cst4dot=diff(cst4)
cst5dot=diff(cst5)
cst6dot=diff(cst6)


syms t Ixx Iyy Izz Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 omega mp mck mst
%% inertia constant of moving platform,cranks & struts
Ip=[Ixx 0 0;0 Iyy 0;0 0 Izz]
Ic=[Ixx1 0 0;0 Iyy1 0;0 0 Izz1]
Ist=[Ixx2 0 0;0 Iyy2 0;0 0 Izz2]

%% Kinetic energy for moving platform
A=transpose(Rxyz)*dRxyz;
Wp=[D(3,2);D(1,3);D(2,1)]

P=[x;y;z]
dotP=[diff(x);diff(y);diff(z)]

KEp=0.5*mp*transpose(dotP)*dotP + 0.5*transpose(Wp)*Ip*Wp
simplify(KEp)



%% kinetic energy For struts %%
R1=Rz1*Ry1*Rx1
omega1=transpose(R1)*diff(R1)
omega1a=omega1(t)
Wst1=[omega1a(3,2);omega1a(1,3);omega1a(2,1)]
KEst1=simplify(0.5*mst*sum(cst1dot.^2) + 0.5*transpose(Wst1)*Ist*Wst1)
simplify(KEst1)

R2=Rz2*Ry2*Rx2
omega2=transpose(R2)*diff(R2)
omega2a=omega2(t)
Wst2=[omega2a(3,2);omega2a(1,3);omega2a(2,1)]
KEst2=simplify(0.5*mst*sum(cst2dot.^2) + 0.5*transpose(Wst2)*Ist*Wst2)
simplify(KEst2)

R3=Rz3*Ry3*Rx3
omega3=transpose(R3)*diff(R3)
omega3a=omega3(t)
Wst3=[omega3a(3,2);omega3a(1,3);omega3a(2,1)]
KEst3=simplify(0.5*mst*sum(cst3dot.^2) + 0.5*transpose(Wst3)*Ist*Wst3)
simplify(KEst3)


R4=Rz4*Ry4*Rx4
omega4=transpose(R4)*diff(R4)
omega4a=omega4(t)
Wst4=[omega4a(3,2);omega4a(1,3);omega4a(2,1)]
KEst4=simplify(0.5*mst*sum(cst4dot.^2) + 0.5*transpose(Wst4)*Ist*Wst4)
simplify(KEst4)

R5=Rz5*Ry5*Rx5
omega5=transpose(R5)*diff(R5)
omega5a=omega5(t)
Wst5=[omega5a(3,2);omega5a(1,3);omega5a(2,1)]
KEst5=simplify(0.5*mst*sum(cst5dot.^2) + 0.5*transpose(Wst5)*Ist*Wst5)
simplify(KEst5)

R6=Rz6*Ry6*Rx6
omega6=transpose(R6)*diff(R6)
omega6a=omega6(t)
Wst6=[omega6a(3,2);omega6a(1,3);omega6a(2,1)]
KEst6=simplify(0.5*mst*sum(cst6dot.^2) + 0.5*transpose(Wst6)*Ist*Wst6)
simplify(KEst6)

%% Total kinetic energy of the system
T=(KEst1+KEst2+KEst3+KEst4+KEst5+KEst6)+(KEck1+KEck2+KEck3+KEck4+KEck5+KEck6)+KEp
T=subs(T,[Vap,Vbp,Vcp],[diff(a(t), t),diff(b(t), t),diff(c(t), t)])

%% constraints of the system from lagrange equation
dtdLdottheta1=diff(diff(T,diff(th1)),t)
dtdLdottheta2=diff(diff(T,diff(th2)),t)
dtdLdottheta3=diff(diff(T,diff(th3)),t)
dtdLdottheta4=diff(diff(T,diff(th4)),t)
dtdLdottheta5=diff(diff(T,diff(th5)),t)
dtdLdottheta6=diff(diff(T,diff(th6)),t)

dtdLdotphi1=diff(diff(T,diff(phi1)),t)
dtdLdotphi2=diff(diff(T,diff(phi2)),t)
dtdLdotphi3=diff(diff(T,diff(phi3)),t)
dtdLdotphi4=diff(diff(T,diff(phi4)),t)
dtdLdotphi5=diff(diff(T,diff(phi5)),t)
dtdLdotphi6=diff(diff(T,diff(phi6)),t)

dtdLdotshi1=diff(diff(T,diff(shi1)),t)
dtdLdotshi2=diff(diff(T,diff(shi2)),t)
dtdLdotshi3=diff(diff(T,diff(shi3)),t)
dtdLdotshi4=diff(diff(T,diff(shi4)),t)
dtdLdotshi5=diff(diff(T,diff(shi5)),t)
dtdLdotshi6=diff(diff(T,diff(shi6)),t)

dtdLdotx=diff(diff(T,diff(x)),t)
dtdLdoty=diff(diff(T,diff(y)),t)
dtdLdotz=diff(diff(T,diff(z)),t)
dtdLdota=diff(diff(T,diff(a)),t)
dtdLdotb=diff(diff(T,diff(b)),t)
dtdLdotc=diff(diff(T,diff(c)),t)

syms x y z ath1 ath2 ath3 ath4 ath5 ath6 aphi1 aphi2 aphi3 aphi4 aphi5 aphi6 ashi1 ashi2 ashi3 ashi4 ashi5 ashi6 ax ay az aa ab ac
vec1=[diff(th1, t, t),diff(th2, t, t),diff(th3, t, t),diff(th4, t, t),diff(th5, t, t),diff(th6, t, t),diff(phi1, t, t),diff(phi2, t, t),diff(phi3, t, t),diff(phi4, t, t),diff(phi5, t, t),diff(phi6, t, t),diff(shi1, t, t),diff(shi2, t, t),diff(shi3, t, t),diff(shi4, t, t),diff(shi5, t, t),diff(shi6, t, t),diff(x, t, t),diff(y, t, t),diff(z, t, t),diff(a, t, t),diff(b, t, t),diff(c, t, t)]
vec2=[ath1,ath2,ath3,ath4,ath5,ath6,aphi1,aphi2,aphi3,aphi4,aphi5,aphi6,ashi1,ashi2,ashi3,ashi4,ashi5,ashi6,ax,ay,az,aa,ab,ac]
zerovec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
%% coefficient of ath1
temp1=subs(expand(dtdLdottheta1),vec1,vec2)
temp2=subs(temp1,vec2,zerovec)
temp3=coeffs(simplify(temp1-temp2),vec2)

%% Mass matrices
[M1,cf1]=coeffs(simplify(subs(expand(dtdLdottheta1),vec1,vec2)-subs(subs(expand(dtdLdottheta1),vec1,vec2),vec2,zerovec)),vec2)
[M2,cf2]=coeffs(simplify(subs(expand(dtdLdottheta2),vec1,vec2)-subs(subs(expand(dtdLdottheta2),vec1,vec2),vec2,zerovec)),vec2)
[M3,cf3]=coeffs(simplify(subs(expand(dtdLdottheta3),vec1,vec2)-subs(subs(expand(dtdLdottheta3),vec1,vec2),vec2,zerovec)),vec2)
[M4,cf4]=coeffs(simplify(subs(expand(dtdLdottheta4),vec1,vec2)-subs(subs(expand(dtdLdottheta4),vec1,vec2),vec2,zerovec)),vec2)
[M5,cf5]=coeffs(simplify(subs(expand(dtdLdottheta5),vec1,vec2)-subs(subs(expand(dtdLdottheta5),vec1,vec2),vec2,zerovec)),vec2)
[M6,cf6]=coeffs(simplify(subs(expand(dtdLdottheta6),vec1,vec2)-subs(subs(expand(dtdLdottheta6),vec1,vec2),vec2,zerovec)),vec2)

[M7,cf7]=coeffs(simplify(subs(expand(dtdLdotphi1),vec1,vec2)-subs(subs(expand(dtdLdotphi1),vec1,vec2),vec2,zerovec)),vec2)
[M8,cf8]=coeffs(simplify(subs(expand(dtdLdotphi2),vec1,vec2)-subs(subs(expand(dtdLdotphi2),vec1,vec2),vec2,zerovec)),vec2)
[M9,cf9]=coeffs(simplify(subs(expand(dtdLdotphi3),vec1,vec2)-subs(subs(expand(dtdLdotphi3),vec1,vec2),vec2,zerovec)),vec2)
[M10,cf10]=coeffs(simplify(subs(expand(dtdLdotphi4),vec1,vec2)-subs(subs(expand(dtdLdotphi4),vec1,vec2),vec2,zerovec)),vec2)
[M11,cf11]=coeffs(simplify(subs(expand(dtdLdotphi5),vec1,vec2)-subs(subs(expand(dtdLdotphi5),vec1,vec2),vec2,zerovec)),vec2)
[M12,cf12]=coeffs(simplify(subs(expand(dtdLdotphi6),vec1,vec2)-subs(subs(expand(dtdLdotphi6),vec1,vec2),vec2,zerovec)),vec2)
  
[M13,cf13]=coeffs(simplify(subs(expand(dtdLdotshi1),vec1,vec2)-subs(subs(expand(dtdLdotshi1),vec1,vec2),vec2,zerovec)),vec2)
[M14,cf14]=coeffs(simplify(subs(expand(dtdLdotshi2),vec1,vec2)-subs(subs(expand(dtdLdotshi2),vec1,vec2),vec2,zerovec)),vec2)
[M15,cf15]=coeffs(simplify(subs(expand(dtdLdotshi3),vec1,vec2)-subs(subs(expand(dtdLdotshi3),vec1,vec2),vec2,zerovec)),vec2)
[M16,cf16]=coeffs(simplify(subs(expand(dtdLdotshi4),vec1,vec2)-subs(subs(expand(dtdLdotshi4),vec1,vec2),vec2,zerovec)),vec2)
[M17,cf17]=coeffs(simplify(subs(expand(dtdLdotshi5),vec1,vec2)-subs(subs(expand(dtdLdotshi5),vec1,vec2),vec2,zerovec)),vec2)
[M18,cf18]=coeffs(simplify(subs(expand(dtdLdotshi6),vec1,vec2)-subs(subs(expand(dtdLdotshi6),vec1,vec2),vec2,zerovec)),vec2)

[M19,cf19]=coeffs(simplify(subs(expand(dtdLdotx),vec1,vec2)-subs(subs(expand(dtdLdotx),vec1,vec2),vec2,zerovec)),vec2)
[M20,cf20]=coeffs(simplify(subs(expand(dtdLdoty),vec1,vec2)-subs(subs(expand(dtdLdoty),vec1,vec2),vec2,zerovec)),vec2)
[M21,cf21]=coeffs(simplify(subs(expand(dtdLdotz),vec1,vec2)-subs(subs(expand(dtdLdotz),vec1,vec2),vec2,zerovec)),vec2)
[M22,cf22]=coeffs(simplify(subs(expand(dtdLdota),vec1,vec2)-subs(subs(expand(dtdLdota),vec1,vec2),vec2,zerovec)),vec2)
[M23,cf23]=coeffs(simplify(subs(expand(dtdLdotb),vec1,vec2)-subs(subs(expand(dtdLdotb),vec1,vec2),vec2,zerovec)),vec2)
[M24,cf24]=coeffs(simplify(subs(expand(dtdLdotc),vec1,vec2)-subs(subs(expand(dtdLdotc),vec1,vec2),vec2,zerovec)),vec2)


cst7=ck7+Ry1*Rx1*[0;0;lst]
cst8=ck8+Ry2*Rx2*[0;0;lst]
cst9=ck9+Ry3*Rx3*[0;0;lst]
cst10=ck10+Ry4*Rx4*[0;0;lst]
cst11=ck11+Ry5*Rx5*[0;0;lst]
cst12=ck12+Ry6*Rx6*[0;0;lst]


ceq1=P+Rxyz*t1-cst7
ceq2=P+Rxyz*t2-cst8
ceq3=P+Rxyz*t3-cst9
ceq4=P+Rxyz*t4-cst10
ceq5=P+Rxyz*t5-cst11
ceq6=P+Rxyz*t6-cst12

syms th1a th2a th3a th4a th5a th6a phi1a phi2a phi3a phi4a phi5a phi6a shi1a shi2a shi3a shi4a shi5a shi6a xa ya za aa ba capos

ceq=[ceq1;ceq2;ceq3;ceq4;ceq5;ceq6];
pos=[th1 th2 th3 th4 th5 th6 phi1 phi2 phi3 phi4 phi5 phi6 shi1 shi2 shi3 shi4 shi5 shi6 x y z a b c];

M1t=M1(t);
M2t=M2(t);
M3t=M3(t);
M4t=M4(t);
M5t=M5(t);
M6t=M6(t);
M7t=M7(t);
M8t=M8(t);
M9t=M9(t);
M10t=M10(t);
M11t=M11(t);
M12t=M12(t);
M13t=M13(t);
M14t=M14(t);
M15t=M15(t);
M16t=M16(t);
M17t=M17(t);
M18t=M18(t);
M19t=M19(t);
M20t=M20(t);
M21t=M21(t);
M22t=M22(t);
M23t=M23(t);
M24t=M24(t);

diagth=[M1t(1);M2t(1);M3t(1);M4t(1);M5t(1);M6t(1)];
thmat=diag(diagth)

diagphi=[M1t(2);M2t(2);M3t(2);M4t(2);M5t(2);M6t(2)];
phimat=diag(diagphi)

diagshi=[M1t(3);M2t(3);M3t(3);M4t(3);M5t(3);M6t(3)];
shimat=diag(diagshi)

diagth1=[M7t(1);M8t(1);M9t(1);M10t(1);M11t(1);M12t(1)]
thmat1=diag(diagth1)

diagphi1=[M7t(2);M8t(2);M9t(2);M10t(2);M11t(2);M12t(2)];
phimat1=diag(diagphi1)


diagth1a=[M13t(1);M14t(1);M15t(1);M16t(1);M17t(1);M18t(1)];
thmat1a=diag(diagth1a)



diagshi1a=[M13t(2);M14t(2);M15t(2);M16t(2);M17t(2);M18t(2)];
shimat1a=diag(diagshi1a)



zeromat=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0];

W=[thmat phimat shimat zeromat];
Q=[thmat1 phimat1 zeromat zeromat];
S=[thmat1a zeromat shimat1a zeromat];

WQS=[W;Q;S];
xvec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mp 0 0 0 0 0];
yvec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mp 0 0 0 0];
zvec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mp 0 0 0];
avec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 M22t(1) M22t(2) M22t(3)];
bvec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 M23t(1) M23t(2) 0];
cvec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 M24t(1) 0 M24t(2)];

M=[WQS;xvec;yvec;zvec;avec;bvec;cvec]
jac = jacobian(ceq,pos)
jacdot=diff(jac);
Jdot1=subs(jacdot,[diff(th1),diff(th2),diff(th3),diff(th4),diff(th5),diff(th6)],[vth1,vth2,vth3,vth4,vth5,vth6])
Jdot2=subs(Jdot1,[diff(phi1),diff(phi2),diff(phi3),diff(phi4),diff(phi5),diff(phi6)],[vphi1,vphi2,vphi3,vphi4,vphi5,vphi6])
Jdot3=subs(Jdot2,[diff(shi1),diff(shi2),diff(shi3),diff(shi4),diff(shi5),diff(shi6)],[vshi1,vshi2,vshi3,vshi4,vshi5,vshi6])
Jdot4=subs(Jdot3,[diff(x),diff(y),diff(z),diff(a),diff(b),diff(c)],[vx,vy,vz,va,vb,vc])

%dq1=[diff(th1); diff(th2); diff(th3); diff(th4); diff(th5); diff(th6); diff(phi1); diff(phi2); diff(phi3); diff(phi4); diff(phi5); diff(phi6); diff(shi1); diff(shi2); diff(shi3); diff(shi4); diff(shi5); diff(shi6); diff(x); diff(y); diff(z); diff(a); diff(b); diff(c)]



%%coriolis matrix
% Define the vector of generalized positions q
q = [th1; th2; th3; th4; th5; th6; phi1; phi2; phi3; phi4; phi5; phi6; shi1; shi2; shi3; shi4; shi5; shi6; x; y; z; a; b; c];
dq=[vth1;vth2;vth3;vth4;vth5;vth6;vphi1;vphi2;vphi3;vphi4;vphi5;vphi6;vshi1;vshi2;vshi3;vshi4;vshi5;vshi6;vx;vy;vz;va;vb;vc];
dM = jacobian(M(:), q); % Compute the Jacobian of the mass matrix
dM1=dM(t);
dMdq = sym(zeros(size(M))); % Initialize Coriolis matrix
for i = 1:24
    for j = 1:24
        for k = 1:24

            dMdq(i, j) = dMdq(i, j) + 0.5 * (dM1(i, j) + dM1(i, k) - dM1(k, j)) * dq(k);
        end
    end
end
C=dMdq

Jdot1=subs(jacdot,[diff(th1),diff(th2),diff(th3),diff(th4),diff(th5),diff(th6)],[vth1,vth2,vth3,vth4,vth5,vth6]);
Jdot2=subs(Jdot1,[diff(phi1),diff(phi2),diff(phi3),diff(phi4),diff(phi5),diff(phi6)],[vphi1,vphi2,vphi3,vphi4,vphi5,vphi6]);
Jdot3=subs(Jdot2,[diff(shi1),diff(shi2),diff(shi3),diff(shi4),diff(shi5),diff(shi6)],[vshi1,vshi2,vshi3,vshi4,vshi5,vshi6]);
Jdot4=subs(Jdot3,[diff(x),diff(y),diff(z),diff(a),diff(b),diff(c)],[vx,vy,vz,va,vb,vc]);
%% Total potential energy

ck=[ck1 ck2 ck3 ck4 ck5 ck6]
cst=[cst1 cst2 cst3 cst4 cst5 cst6]


cst1=ck7+Ry1*Rx1*[0;0;lst/2]
ck1=Rz1*[rb;(lc/2)*cos(th1);(lc/2)*sin(th1)]

zck1=ck1(t);
zcst1=cst1(t);
zck2=ck2(t);
zcst2=cst2(t);
zck3=ck3(t);
zcst3=cst3(t);
zcst4=cst4(t);
zck4=ck4(t);
zck5=ck5(t);
zcst5=cst5(t);
zck6=ck6(t);
zcst6=cst6(t);

g=9.8;
PE1= mck*g*zck1(3)+mst*g*zcst1(3); 
PE2= mck*g*zck2(3)+mst*g*zcst2(3); 
PE3= mck*g*zck3(3)+mst*g*zcst3(3); 
PE4= mck*g*zck4(3)+mst*g*zcst4(3); 
PE5= mck*g*zck5(3)+mst*g*zcst5(3); 
PE6= mck*g*zck6(3)+mst*g*zcst6(3); 


PE=PE1+PE2+PE3+PE4+PE5+PE6+mp*g*z;
%% total energy of the system
TE=T+PE;
TE1=subs(TE,[diff(th1),diff(th2),diff(th3),diff(th4),diff(th5),diff(th6)],[vth1,vth2,vth3,vth4,vth5,vth6]);
TE2=subs(TE1,[diff(phi1),diff(phi2),diff(phi3),diff(phi4),diff(phi5),diff(phi6)],[vphi1,vphi2,vphi3,vphi4,vphi5,vphi6]);
TE3=subs(TE2,[diff(shi1),diff(shi2),diff(shi3),diff(shi4),diff(shi5),diff(shi6)],[vshi1,vshi2,vshi3,vshi4,vshi5,vshi6]);
TE4=subs(TE3,[diff(x),diff(y),diff(z),diff(a),diff(b),diff(c)],[vx,vy,vz,va,vb,vc])


% Assuming you have defined all your symbolic variables and expressions as in your original code




% G=jacobian(PE(:),q);
% Minv = M \ eye(size(M));
%  A_prime= jac*Minv*transpose(jac);
% % %%Ainv=A_prime \ eye(size(A_prime))
% % F=-C*dq-transpose(G)
% % %%Lmda=-Ainv*(jacdot*dq+jac*Minv*F)
% % E=jac*Minv*F;
% % S=jacdot*dq;
% % %Ainv=(1/det(A_prime))*adjugate(A_prime);
% % [L,U,P]=lu(A_prime);
% % A_inv=inv(U)*inv(L)*P
% 
% fid = fopen('mass.txt', 'wt');
% fprintf(fid, '%s\n', char(M));
% fclose(fid);
% 
% fid = fopen('coriolis.txt', 'wt');
% fprintf(fid, '%s\n', char(C));
% fclose(fid);
% 
% fid = fopen('jac.txt', 'wt');
% fprintf(fid, '%s\n', char(jac));
% fclose(fid);
% 
% fid = fopen('jacdot.txt', 'wt');
% fprintf(fid, '%s\n', char(Jdot4));
% fclose(fid);
% 
% fid = fopen('Gravitatinal.txt', 'wt');
% fprintf(fid, '%s\n', char(G));
% fclose(fid);
% 
% fid = fopen('Total_energy.txt', 'wt');
% fprintf(fid, '%s\n', char(TE4));
% fclose(fid);


%writematrix(M)