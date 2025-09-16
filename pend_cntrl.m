function dydt = pend_cntrl(t,y)
%simpend control
%
des=pi;
desv=0;
g=9.8;
L=1;
kp=300;
kv=5;

dydt = [y(2); -(g/L)*y(1)+kp*(des-y(1))+kv*(desv-y(2))];