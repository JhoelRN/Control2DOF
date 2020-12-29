%% SCA - GRUPO A
%% CONTROL 2 DOF BRAZO ROBOTICO
%% PASSIVITY BASED CONTROL ADAPTATIV
clc; clear all; close all

%Parametros del brazo
global m1 l1 l2 m2t g pi
m1 =1; l1 =1; l2=1; g = 9.8; pi = 3.14;
%Tiempos
tf = 40; 
global param;

global ts = [0 tf];

global a;  %a es tiempo aux
global torque error mass
torque =[];
error=[];
mass =[];
%condiciones iniciales
%x0 = [0.05,0.1,0.05,0.1,5]; %principal
x0 = [1.5,0.1,0.05,0.1,30];%extraCases
% cond iniciales : theta1, theta2, . . . , masa

function [dx] = planarArmODEadaptive(t,x)
global m2t mass
param = 1.8;

%Case 1 - (P-0.2, lambda = 0.99 tf =30)
%m2t = sin(t) + 3; %(SENOIDAL)

%Case 2 - (P-1.5, lambda= 0.99 tf = 20)
%m2t = 2; %OK!  (ESTABILIZACION CONSTANTE)

%Case 3 - (P-1.5, lambda= 0.99 tf = 20) 
% if t> 4   / param     % (ESTABILIZACION ESCALON)
%     m2t = 2;
%     mass =[mass, m2t];
% elseif t<8  / param
%     m2t = 4;
%     mass = [mass, m2t];
% end
 
 
u=(t>=0); 
%Case 4 -  lambda =0.99 tf = 20
if t> 0.0001     % (ESTABILIZACION RAMPA)
     m2t = t.*u;
     mass =[mass, m2t];
 elseif t<0.0001
     m2t = 0.01;
     mass = [mass, m2t];
end
 
% trayectorias deseadas
theta_d = [0;sin(t)];
dtheta_d = [0;cos(t)];%primera derivada
ddtheta_d = [0;-sin(t)];%segunda derivada

%  trayectorias
theta = x(1:2);
dtheta= x(3:4);

% errores
global lambda e de a v r m1 l1 l2 g
lambda = 0.99;
e = theta - theta_d;
de = dtheta - dtheta_d;
a = ddtheta_d - (lambda*de);
v = dtheta_d - (lambda*e);
r = de + (lambda*e);

%error tracking  LISTA..
global error
error = [error, e];

% Matriz definida positiva (se utilizará para W_update)
P = 1.5*eye(11);

% Modelo
global M1 G1 M2 C2 G2 M C G PM PG PC
% El modelo dinámico del sistema 
%-.- se caracteriza por M y C

%para link 1
M1 = [(1/3)*m1*(l1)^2 , 0 ; 0, 0];
G1 = [(1/2)*m1*g*l1*cos(x(1)); 0];

%para link 2
PM = [l1^2 + (1/3)*l2^2 + l1*l2*cos(x(2)) , (1/3)*l2^2 + (1/2)*l1*l2*cos(x(2)); (1/3)*l2^2 + (1/2)*l1*l2*cos(x(2)), (1/3)*l2^2];
PC = [0, -((1/2)*l1*l2*sin(x(2))*x(4) + l1*l2*sin(x(2))*x(3)) ; (1/2)*l1*l2*sin(x(2))*x(3), 0];
PG = [(1/2)*g*l1*cos(x(1)) + (1/2)*g*l2*cos(x(1)+ x(2)) ; (1/2)*g*l2*cos(x(1)+ x(2))];
M2 = m2t*PM;
G2 = m2t*PG;
C2 = m2t*PC;

%modelo actual
M = M1 + M2;
C = C2;
G = G1 + G2;
invM = inv(M);
invMC = inv(M)*C;

% Series Fourier 
Z = [(1/2); cos((pi*(t))/5);sin((pi*(t))/5);cos((2*pi*(t))/5);sin((2*pi*(t))/5);cos((3*pi*(t))/5);sin((3*pi*(t))/5);cos((4*pi*(t))/5);sin((4*pi*(t))/5);cos((5*pi*(t))/5);sin((5*pi*(t))/5)];

m2t_bar = x(5);

% Modelo estimado
global M_bar C_bar G_bar M2_bar C2_bar G2_bar

M2_bar = m2t_bar*[l1^2 + (1/3)*l2^2 + l1*l2*cos(x(2)) , (1/3)*l2^2 + (1/2)*l1*l2*cos(x(2));
                    (1/3)*l2^2 + (1/2)*l1*l2*cos(x(2)), (1/3)*l2^2];
                
C2_bar = m2t_bar*[0,                -((1/2)*l1*l2*sin(x(2))*x(4) + l1*l2*sin(x(2))*x(3)) ;
         (1/2)*l1*l2*sin(x(2))*x(3),                                    0];
     
G2_bar = m2t_bar*[(1/2)*g*l1*cos(x(1)) + (1/2)*g*l2*cos(x(1)+ x(2)) ; (1/2)*g*l2*cos(x(1)+ x(2))];

M_bar = M1 + M2_bar;
C_bar = C2_bar;
G_bar = G1 + G2_bar;

% Torque
tau = adaptive_ctrl(theta_d, dtheta_d, ddtheta_d, theta, dtheta);
global torque
torque = [torque, tau];

%Actualizacion estdo del sistema, 
%calculo de dx
dx=zeros(5,1);
dx(1) = x(3);
dx(2) = x(4);
dx(3:4) = -invMC* x(3:4) - invM*G + invM*tau; % porque ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
%calculo ley - W_update (actualizacion)
W_update = -inv(P)*[Z*transpose(e)*PM*a + Z*transpose(e)*PC*v + Z*transpose(e)*PG];%% ec(6)
dx(5) = transpose(W_update)*Z;

end

% Calculo de Torque
function tau = adaptive_ctrl(theta_d, dtheta_d, ddtheta_d, theta, dtheta)
global M C M_bar C_bar G_bar lambda e de a v r
%Kp = 100*eye(1);
Kv = 500*eye(2);
tau = (M_bar*a)+ (C_bar*v) + (G_bar) - Kv*r;
end


% Implementacion Control adaptativo pasiv
%options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);
[T,X] = ode23(@(t,x)planarArmODEadaptive(t,x),ts,x0);


% %Case 4 -  lambda =0.99 tf = 20
u=(a>=0); 
a = 0:0.5:40;
if a> 0.001    % (ESTABILIZACION RAMPA)
     m2t = a.*u;
elseif a<0.001
     m2t = 0.01;
end

figure('name',' Theta-1 Control Adaptativo');
plot(T, X(:,1),'r-');
hold on
plot(T, 0*ones(size(T,1),1), 'b-');
legend('Theta Actual', 'theta deseado')
xlabel('Tiempo')
ylabel('Theta 1')
title('Theta-1 Control Adaptativo');
print -djpg fig4_1.jpg %para guardar la imagen (formato.jpg)


figure('name','Theta-2 Control Adaptativo');
plot(T, X(:,2),'r-');
hold on
plot(T, sin(T), 'b-');
legend('Theta Actual', 'theta deseado')
xlabel('Tiempo')
ylabel('Theta 2')
title('Theta-2 Control adaptativo');
print -djpg fig4_2.jpg %para guardar la imagen (formato.jpg)

figure('name','Error Masa');
plot(T, X(:,5),'r-');
hold on
%plot(T, mass(1, 1: size(T,1)) ,'b-');%caso escalon
%
%plot(T, 3+sin(T),'b-');%caso sen
%
%plot(T, 2*ones(size(T,1),1),'b-');%caso cte
%
plot(a, a ,'m-');%caso rampa
xlabel('Tiempo');
ylabel('Masa');
legend('Masa estimada', 'Masa actual:2');
title('Error masa - Control Adaptativo');
print -djpg fig4_3.jpg %para guardar la imagen (formato.jpg)




%plot(T, X(:,5),'r-');
%plot(a,a,'r-');
