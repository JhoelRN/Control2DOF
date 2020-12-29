% ---------- Diseño del sistema regulador óptimo cuadrático ---------
close all; clear;


 A  =[0.          1.         0.          0;      
  -0.009805    0.000012   0.0196836   0.000004;
   0.          0.         0.          1;      
   0.0882748  -0.000068  -0.0984476  -0.00002]; 


 B  = [0.         0;      
   0.999996  -4.999978;
   0.         0;      
  -4.999978   28.99988];
  
 C  = [1.   0.   0.   0;
   0.   1.   0.   0;
   0.   0.   1.   0;
   0.   0.   0.   1.];
    
Q=C'*C;

iden=[1 0; 0 1]; 

rho=0.1;
R=rho*iden;
[G,P,E]=lqr(A,B,Q,R)

%% closed loop systems - Se debde de enteder esta movida pARA tracking
% sys=ss(A-B*G, eye(4),eye(4),eye(4)); 
sys=ss(A-B*G, eye(4),eye(4),eye(4)); 
% impulse response at initial condiction
t=0:0.01:20; 
x=initial(sys,[1;0;1;0],t);

plot(t,x)
legend('th1(rad)','thdot1(rad/seg)','th2(rad)','thdot2(rad/seg)')
grid on;
grid minor;
xlabel('Time (s)')
ylabel('Amplitud de estados') 
title('Respuesta del sistema')








