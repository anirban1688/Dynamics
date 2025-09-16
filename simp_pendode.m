clc
clear all
[t,y] = ode45(@pend_cntrl,[0 50],[1; 5]);
%Plot the solutions for y


plot(t,y(:,1) ,t,y(:,2),'o')
title('Solution of simple pendulum with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('pos',' velocity')

figure
 
L=1;
axlen=2;
x1=L*sin(y(:,1));
y1=-L*cos(y(:,1));
O1=[0 0];
axis(gca,'equal');
axis([-axlen axlen -axlen axlen]); %%%%%%%%%%%% fix the frame dimensions %%%%%%%%
grid on;
    
for i=1:length(t)
    P1=[x1(i,1) y1(i,1)]; %%%%%%%%%%% define the position of the bob%%%%%%%%%%   
    O_circ=viscircles(O1,0.01);%%%%%%%%% draw the origin%%%%%%%%%%%
    %%%%%%%%%%% viscircles requires the image processing toolbox %%%%%%%%%% 
    pend=line([O1(1) P1(1)], [O1(2) P1(2)]);%%%%%% draw the link %%%%%%%%    
    ball =viscircles(P1,0.05); %%%%%%%%% draw the pendulum bob%%%%%%%    
    pause(0.1);
    
    %%%%%%%%%%%%delete the previous objects%%%%%%%%%%%%%%%
    
    if i<length(t)
        delete(pend);
        delete(ball);
        delete(O_circ);
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

