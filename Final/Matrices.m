close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE569

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints / Robot Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix a
syms Theta_a D_a;

a = [cos(Theta_a) -sin(Theta_a)    0              0
     sin(Theta_a)  cos(Theta_a)    0              0
     0              0              1              D_a
     0              0              0              1];


% Matrix b
syms D_b

b = [1              0              0              0
     0              0              1              0
     0             -1              0              D_b              
     0              0              0              1];


% Matrix f
syms Theta_f A_f;

f = [cos(Theta_f) -sin(Theta_f)    0              A_f*cos(Theta_f)
     sin(Theta_f)  cos(Theta_f)    0              A_f*sin(Theta_f)
     0              0              1              0              
     0              0              0              1];


% Analitically, look at the results of Tabd (or T3)
% We found that sin(Theta_f) = 0 and cos(Theta_f)=1
% This is in line if we assume that no changes in alpha took place
%      between Matrices b and f

% This is noted on Manual Calculations done in the test
f = [1              0              0              A_f
     0              1              0              0
     0              0              1              0              
     0              0              0              1];


% Matrix for Frankenbot
%Tabf = a*b*f
Tabf = a*b*f;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is based on Manual Calculations done in the test

Jv1 = [A_f*sin(Theta_a)   0                  0
      -A_f*cos(Theta_a)   0                  0
       0                  0                  0];

Jv1i = [A_f*sin(Theta_a) -A_f*cos(Theta_a)   0
        0                 0                  0
        0                 0                  0];


Jv2 = [A_f*sin(Theta_a)   A_f*sin(Theta_a)   0
      -A_f*cos(Theta_a)  -A_f*cos(Theta_a)   0
       0                  0                  0];


Jv2i = [A_f*sin(Theta_a) -A_f*cos(Theta_a)   0
        A_f*sin(Theta_a) -A_f*cos(Theta_a)   0
        0                 0                  0];


Jv3 = [A_f*sin(Theta_a)   A_f*sin(Theta_a)  -sin(Theta_a)
      -A_f*cos(Theta_a)  -A_f*cos(Theta_a)   cos(Theta_a)
       0                  0                  0];

Jv3i = [A_f*sin(Theta_a) -A_f*cos(Theta_a)   0
        A_f*sin(Theta_a) -A_f*cos(Theta_a)   0
       -sin(Theta_a)      cos(Theta_a)       0];


Jw1 = [0  0  0
       0  0  0
       1  0  0];

Jw1i = [0  0  1
        0  0  0
        0  0  0];


Jw2 = [0  0  0
       0  0  0
       1  1  0];

Jw2i = [0  0  1
        0  0  1
        0  0  0];


Jw3 = [0  0  0
       0  0  0
       1  1  0];

Jw3i = [0  0  1
        0  0  1
        0  0  0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is based on Manual Calculations done in the test

R1 = [cos(Theta_a) -sin(Theta_a) 0
      sin(Theta_a)  cos(Theta_a) 0
      0             0            1];

R1i = [cos(Theta_a)  sin(Theta_a) 0
     -sin(Theta_a)  cos(Theta_a) 0
      0             0            1];

R2 = [cos(Theta_a)  0 -sin(Theta_a)
      sin(Theta_a)  0  cos(Theta_a) 
      0            -1  0];

R2i = [cos(Theta_a)  sin(Theta_a)  0
       0             0            -1 
      -sin(Theta_a)  cos(Theta_a)  0];

R3 = [cos(Theta_a)  0 -sin(Theta_a)
      sin(Theta_a)  0  cos(Theta_a) 
      0            -1  0];

R3i = [cos(Theta_a)  sin(Theta_a)  0
       0             0            -1 
      -sin(Theta_a)  cos(Theta_a)  0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertia Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is based on Manual Calculations done in the test

syms m1 m2 m3 I1 I2 I3

 
I1  = [(5/12)*m1*(D_a^2)  0                  0
       0                  (5/12)*m1*(D_a^2)  0
       0                  0                  (1/6)*m1*(D_a^2)];

I2  = [(5/12)*m2*(D_b^2)  0                  0
       0                  (5/12)*m2*(D_b^2)  0
       0                  0                  (1/6)*m2*(D_b^2)];

I3  = [(1/6)*m3*(A_f^2)  0                  0
       0                 (5/12)*m3*(A_f^2)  0
       0                 0                  (5/12)*m3*(A_f^2)];


D3 = m1*(Jv1i)*Jv1 + (Jw1i)*R1*I1*(R1i)*Jw1 ...
   + m2*(Jv2i)*Jv2 + (Jw2i)*R2*I2*(R2i)*Jw2 ...
   + m3*(Jv3i)*Jv3 + (Jw3i)*R3*I3*(R3i)*Jw3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining Generalized Coordonates & Derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms Theta_a_dot A_f_dot Theta_a_double_dot A_f_double_dot

% Generalized Coordonates
q  = [Theta_a; 3*pi/2; A_f];
q_dot = [Theta_a_dot; 0 ;A_f_dot];
q_double_dot = [Theta_a_double_dot; 0 ;A_f_double_dot];

% Transposes
qi = [Theta_a; 3*pi/2; A_f];
qi_dot = [Theta_a_dot 0 A_f_dot];
qi_double_dot = [Theta_a_double_dot 0 A_f_double_dot];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the Kinematic Energy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = (1/2)*qi_dot*D3*q_dot;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the C(q,q_ador) Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is based on Manual Calculations done in the test

% Splitting in columns for ease of readability
C_col1=[(A_f*A_f_dot*(m1+m2+(17/12)*m3))
        (A_f*A_f_dot*(m2+(17/12)*m3))
        A_f*Theta_a_dot*(m1+m2+(17/12)*m3)-A_f*A_f_dot*m3];

C_col2 = [(A_f*A_f_dot*(m2+(17/12)*m3))
          (A_f*A_f_dot*(m2+(17/12)*m3))
          (A_f*Theta_a_dot*(m2+(17/12)*m3)-A_f*A_f_dot*m3)];

C_col3 = [-(A_f*Theta_a_dot*(m1+m2+(17/12)*m3))
           (A_f*Theta_a_dot*(m2+(17/12)*m3))
           ((-1/2)*Theta_a_dot*m3+A_f*Theta_a_dot*m3)];

C = [C_col1 C_col2 C_col3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the Gravity Vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is based on Manual Calculations done in the test
% Note: This makes since,a s the robot does not move on the z-axis
g_p = zeros(3,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the generalized forces:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: D(q)*q_double_dot + C(q,q_dot)*q_dot + g(p) = Tau (Eq. 6.66)

syms Tau
Tau = D3*q_double_dot + C*q_dot + g_p;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paramter required to graph Tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Robot Mesurements in inches
D_a_num = 14;
D_b_num = 18;

% Update Tau with values
Tau_num = subs(Tau,D_a,D_a_num);
Tau_num = subs(Tau_num,D_b,D_b_num);

rho = 0.1; % Aluminium density in lb / in^3

m1_num = 2*rho*(D_a_num^3);
m2_num = 2*rho*(D_b_num^3);
m3_num = 2*rho*(A_f^3);

% Update Tau with values
Tau_num = subs(Tau_num,m1,m1_num);
Tau_num = subs(Tau_num,m2,m2_num);
Tau_num = subs(Tau_num,m3,m3_num); % Update for A_f will be done with q

% Generalized Coordonates

% We will assume Theta_a follows the a behaviour of the shape
%       y = -(x^2)+2*pi*x
syms Theta_a_num
Theta_a_num = -(Theta_a^2)+2*pi*Theta_a;
Theta_a_dot_num = diff(Theta_a_num, Theta_a);
Theta_a_double_dot_num = diff(Theta_a_dot_num, Theta_a);

% Update Tau with values
Tau_num = subs(Tau_num,Theta_a,Theta_a_num);
Tau_num = subs(Tau_num,Theta_a_dot,Theta_a_dot_num);
Tau_num = subs(Tau_num,Theta_a_double_dot,Theta_a_double_dot_num);


% We will assume A_f follows the a behaviour of the shape
%       y = x + 2
syms A_f_num A_f_dot_num A_f_double_dot_num
A_f_num = A_f + 2;
A_f_dot_num = diff(A_f_num, A_f);
A_f_double_dot_num = diff(A_f_dot_num, Theta_a);

% Update Tau with values
Tau_num = subs(Tau_num,A_f,A_f_num);
Tau_num = subs(Tau_num,A_f_dot,A_f_dot_num);
Tau_num = subs(Tau_num,A_f_double_dot,A_f_double_dot_num);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1 - q(1) - Theta_a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(3,1,1)
plot(rand(1,10),'.-');
title('q(1) - \theta_a') 
subplot(3,1,2)
plot(rand(1, 10),'.-');
title('dq(1) - d\theta_a') 
subplot(3,1,3)
plot(rand(1, 15),'.-');
title('ddq(1) - dd\theta_a') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure X - q(2) - alpha_a 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since alpha is a constant, no graph will be plotted.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3 - q(3) - A_f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
subplot(3,1,1)
plot(rand(1,10),'.-');
title('q(3) - A_f') 
subplot(3,1,2)
plot(rand(1, 10),'.-');
title('dq(3) - dA_f') 
subplot(3,1,3)
plot(rand(1, 15),'.-');
title('ddq(3) - ddA_f') 

