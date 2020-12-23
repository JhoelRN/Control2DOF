%% SCA - GRUPO A
%% CONTROL 2 DOF BRAZO ROBOTICO
%% PASSIVITY BASED ADAPTIVE CONTROL
clc; clear all; close all

%Parametros del brazo
global m1 l1 l2 m2t g pi
m1 =1; l1 =1; l2=1; g = 9.8; pi = 3.14;
%time steps
tf = 20;

global torque error mass
torque =[];
error=[];
mass =[];
%condiciones iniciales
x0 = [0.05,0.1,0.05,0.1,5];


function [dx] = planarArmODEadaptive(t,x)
global m2t mass
% Case 1 - (P-0.2, lambda = 0.99 tf =30)
%m2t = sin(t) + 3;

%Case 2 - (P-1.5, lambda= 0.99 tf = 20)
%m2t = 2; %OK!

%Case 3 - (P-1.5, lambda= 0.99 tf = 20)
 if t< 7
     m2t = 4;
     mass=[mass, m2t];
 elseif t>7
     m2t = 2;
     mass = [mass, m2t];
 end

% trayectorias deseadas
theta_d = [0;sin(t)];
dtheta_d = [0;cos(t)];
ddtheta_d = [0;-sin(t)];

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

%tracking error list
global error
error = [error, e];

% a positive definite matrix (to be used later for W_update)
P = 1.5*eye(11);

% True model
global M1 G1 M2 C2 G2 M C G PM PG PC
% actual dynamic model of the system is characterized by M and C
%for link 1
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

% Fourier Series
Z = [(1/2); cos((pi*(t))/5);sin((pi*(t))/5);cos((2*pi*(t))/5);sin((2*pi*(t))/5);cos((3*pi*(t))/5);sin((3*pi*(t))/5);cos((4*pi*(t))/5);sin((4*pi*(t))/5);cos((5*pi*(t))/5);sin((5*pi*(t))/5)];

m2t_bar = x(5);

% Estimated model
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

%update the system state, compute dx
dx=zeros(5,1);
dx(1) = x(3);
dx(2) = x(4);
dx(3:4) = -invMC* x(3:4) - invM*G + invM*tau; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
%update law - W_update
W_update = -inv(P)*[Z*transpose(e)*PM*a + Z*transpose(e)*PC*v + Z*transpose(e)*PG];%% ec(6)
dx(5) = transpose(W_update)*Z;

end




% function to calculate torque
function tau = adaptive_ctrl(theta_d, dtheta_d, ddtheta_d, theta, dtheta)
global M C M_bar C_bar G_bar lambda e de a v r
%Kp = 100*eye(1);
Kv = 500*eye(2);
tau = (M_bar*a)+ (C_bar*v) + (G_bar) - Kv*r;
end















% Implementacion Control adaptativo pasiv
%options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);
[T,X] = ode45(@(t,x)planarArmODEadaptive(t,x),[0 tf],x0);

figure('name',' Theta-1 Control Adaptativo');
plot(T, X(:,1),'r-');
hold on
plot(T, 0*ones(size(T,1),1), 'b-');
legend('Theta Actual', 'theta deseado')
xlabel('Tiempo')
ylabel('Theta 1')
title('Theta-1 Control Adaptativo');


figure('name','Theta-2 Control Adaptativo');
plot(T, X(:,2),'r-');
hold on
plot(T, sin(T), 'b-');
legend('Theta Actual', 'theta deseado')
xlabel('Tiempo')
ylabel('Theta 2')
title('Theta-2 Control adaptativo');

figure('name','Error Masa');
plot(T, X(:,5),'r-');
hold on
plot(T, mass(1,1:size(T,1)) ,'b-');%caso2
%plot(T, 3+sin(T),'b-');
%plot(T, 2*ones(size(T,1),1),'b-');
xlabel('Tiempo');
ylabel('Masa');
legend('Masa estimada', 'Masa actual:2');
title('Error masa - Control Adaptativo');


