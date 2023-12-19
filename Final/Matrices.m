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

syms rho m1 m2 m3 I1 I2 I3

%rho =

m1_num = 2*rho*(D_a^3);
m2_num = 2*rho*(D_b^3);
m3_num = 2*rho*(A_f^3);
 
I1_num = [(5/12)*m1*(D_a^2)  0                  0
          0                  (5/12)*m1*(D_a^2)  0
          0                  0                  (1/6)*m1*(D_a^2)];

I2_num  = [(5/12)*m2*(D_b^2)  0                  0
           0                  (5/12)*m2*(D_b^2)  0
           0                  0                  (1/6)*m2*(D_b^2)];

I3_num  = [(1/6)*m3*(A_f^2)  0                  0
           0                 (5/12)*m3*(A_f^2)  0
           0                 0                  (5/12)*m3*(A_f^2)];


D3 = m1*(Jv1i)*Jv1 + (Jw1i)*R1*I1*(R1i)*Jw1 ...
   + m2*(Jv2i)*Jv2 + (Jw2i)*R2*I2*(R2i)*Jw2 ...
   + m3*(Jv3i)*Jv3 + (Jw3i)*R3*I3*(R3i)*Jw3;

D3_num = m1*(Jv1i)*Jv1 + (Jw1i)*R1*I1_num*(R1i)*Jw1 ...
       + m2*(Jv2i)*Jv2 + (Jw2i)*R2*I2_num*(R2i)*Jw2 ...
       + m3*(Jv3i)*Jv3 + (Jw3i)*R3*I3_num*(R3i)*Jw3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining Generalized Coordonates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q1 = Theta_a;
q2 = 3*pi/2;
q3 = A_f;

q  = [q1; q2; q3];
qi = [q1 q2 q3];

% Defining Generalized Derivatives Coordonates
syms Theta_a_dot A_f_dot
q_dot = [Theta_a_dot; 0 ;A_f_dot];
qi_dot = [Theta_a_dot 0 A_f_dot];

% Defining Generalized Derivatives 2nd Coordonates
syms Theta_a_double_dot A_f_double_dot
q_double_dot = [Theta_a_double_dot; 0 ;A_f_double_dot];
qi_double_dot = [Theta_a_double_dot 0 A_f_double_dot];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the Kinematic Energy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = (1/2)*qi_dot*D3*q_dot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the C(q,q_ador) Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is based on Manual Calculations done in the test

C_col1 = [(A_f*A_f_dot*(m1_num+m2_num+(17/12)*m3_num))
          (A_f*A_f_dot*(m2_num+(17/12)*m3_num))
          (A_f*Theta_a_dot*(m1_num+m2_num+(17/12)*m3_num)-A_f*A_f_dot*m3_num)]

C_col2 = [(A_f*A_f_dot*(m2_num+(17/12)*m3_num))
          (A_f*A_f_dot*(m2_num+(17/12)*m3_num))
          (A_f*Theta_a_dot*(m2_num+(17/12)*m3_num)-A_f*A_f_dot*m3_num)];

C_col3 = [-(A_f*Theta_a_dot*(m1_num+m2_num+(17/12)*m3_num))
           (A_f*Theta_a_dot*(m2_num+(17/12)*m3_num))
           ((-1/2)*Theta_a_dot*m3_num+A_f*Theta_a_dot*m3_num)];

C = [C_col1 C_col2 C_col3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the Gravity Vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is based on Manual Calculations done in the test


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the generalized forces:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: D(q)*q_double_dot + C(q,q_dot)*q_dot + g(p) = Tau (Eq. 6.66)

syms Tau
Tau = D3*q_double_dot + C













