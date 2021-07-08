%% VERIFY FORWARD KINEMATICS
clear all; clc;

T_01 = SE3(SO3());
T_01.t = [0, 0, 0.07];

theta1 = 0;
T_11dot = SE3(SO3(rotz(theta1)));

T_1dot2 = SE3(SO3([0, 0, 1; 1, 0, 0; 0, 1, 0]));
T_1dot2.t = [0.08, 0, 0.12];

theta2 = 0;
T_22dot = SE3(SO3(rotz(theta2)));

T_2dot3 = SE3(SO3());
T_2dot3.t = [0, 0.43, 0];

theta3 = 0;
T_33dot = SE3(SO3(rotz(theta3)));

T_3dot4 = SE3(SO3([0, -1, 0; 0, 0, 1; -1, 0 ,0]));
T_3dot4.t = [0, 0.335, -0.08];

T_01dot = T_01 * T_11dot;
T_02 = T_01dot * T_1dot2;
T_02dot = T_02 * T_22dot;
T_03 = T_02dot * T_2dot3;
T_03dot = T_03 * T_33dot;
T_04 = T_03dot * T_3dot4;

trplot(T_01, 'length', 0.2, 'color','k')
view(45, 45)
hold on
trplot(T_01dot, 'length', 0.2, 'color','b')
trplot(T_02 , 'length', 0.2, 'color','k')
trplot(T_02dot , 'length', 0.2, 'color','b')
trplot(T_03 , 'length', 0.2, 'color','k')
trplot(T_03dot , 'length', 0.2, 'color','b')
trplot(T_04 , 'length', 0.2, 'color','k')

%% JOINT JACOBIAN
clear all; clc;

syms theta1 theta2 theta3 theta4

T_01 = SE3(SO3());
T_01.t = [0, 0, 0.07];

T_1dot2 = SE3(SO3([0, 0, 1; 1, 0, 0; 0, 1, 0]));
T_1dot2.t = [0.08, 0, 0.12];

T_2dot3 = SE3(SO3());
T_2dot3.t = [0, 0.43, 0];

T_11dot = [cos(theta1) -sin(theta1) 0 0;
           sin(theta1) cos(theta1) 0 0;
           0 0 1 0
           0 0 0 1];

T_22dot = [cos(theta2) -sin(theta2) 0 0;
           sin(theta2) cos(theta2) 0 0;
           0 0 1 0
           0 0 0 1];

T_33dot = [cos(theta3) -sin(theta3) 0 0;
           sin(theta3) cos(theta3) 0 0;
           0 0 1 0
           0 0 0 1];
       

T_01dot = T_01.T * T_11dot;
T_02dot = T_01dot * T_1dot2.T * T_22dot;
T_03dot = T_02dot * T_2dot3.T * T_33dot;

xp = T_03dot(1:3,4);
Jv1 = diff(xp,theta1);
Jv2 = diff(xp,theta2);
Jv3 = diff(xp,theta3);
Jw1 = T_01dot(1:3,3);
Jw2 = T_02dot(1:3,3);
Jw3 = T_03dot(1:3,3);

Jacobian = [Jv1 Jv2 Jv3 ;
            Jw1 Jw2 Jw3 ]
        
%% COM JACOBIAN
clear all; clc;

syms theta1 theta2 theta3 theta4

T_01 = SE3(SO3());
T_01.t = [0, 0, 0.07];

T_1dot2 = SE3(SO3([0, 0, 1; 1, 0, 0; 0, 1, 0]));
T_1dot2.t = [0.08, 0, 0.12];

T_2dot3 = SE3(SO3());
T_2dot3.t = [0, 0.43, 0];


T_11dot = [cos(theta1) -sin(theta1) 0 0;
           sin(theta1) cos(theta1) 0 0;
           0 0 1 0
           0 0 0 1];

T_22dot = [cos(theta2) -sin(theta2) 0 0;
           sin(theta2) cos(theta2) 0 0;
           0 0 1 0
           0 0 0 1];

T_33dot = [cos(theta3) -sin(theta3) 0 0;
           sin(theta3) cos(theta3) 0 0;
           0 0 1 0
           0 0 0 1];
       

T_01dot = T_01.T * T_11dot;
T_02dot = T_01dot * T_1dot2.T * T_22dot;
T_03dot = T_02dot * T_2dot3.T * T_33dot;

T_1dot1com = [1 0 0 0.008355;
              0 1 0 0;
              0 0 1 0.102659;
              0 0 0 1];
T_01com = T_01dot * T_1dot1com;

T_2dot2com = [1 0 0 0;
              0 1 0 0.213502;
              0 0 1 0.0919281;
              0 0 0 1];
T_02com = T_02dot * T_2dot2com;

T_3dot3com = [1 0 0 0;
              0 1 0 0.120975;
              0 0 1 -0.0755858;
              0 0 0 1];
T_03com = T_03dot * T_3dot3com;

          
Pc1 = T_01com(1:3,4);
Jvc1 = [diff(Pc1, theta1) diff(Pc1, theta2) diff(Pc1, theta3)];
          
Pc2 = T_02com(1:3,4);
Jvc2 = [diff(Pc2, theta1) diff(Pc2, theta2) diff(Pc2, theta3)];

Pc3 = T_03com(1:3,4);
Jvc3 = [diff(Pc3, theta1) diff(Pc3, theta2) diff(Pc3, theta3)];
          
Jw1 = [T_01com(1:3,3) zeros(3,1) zeros(3,1)];
Jw2 = [T_01com(1:3,3) T_02com(1:3,3) zeros(3,1)];
Jw3 = [T_01com(1:3,3) T_02com(1:3,3) T_03com(1:3,3)];


m1 = 2.75458;
m2 = 8.143;
m3 = 5.14886;

I1_body = [0.0109697 -2.08287e-08 -0.000398309;
           -2.08287e-08 0.011856 2.03974e-09;
           -0.000398309 2.03974e-09 0.00604953];

       
I1_world = T_01com(1:3,1:3) * I1_body * transpose(T_01com(1:3,1:3));

I2_body = [0.25374 1.41285e-07 -7.13886e-08;
           1.41285e-07 0.0212824 0.00499483;
           -7.13886e-08 0.00499483 0.247523];
       
I2_world = T_02com(1:3,1:3) * I2_body * transpose(T_02com(1:3,1:3));
       
I3_body = [0.075434 2.67415e-08 2.13561e-09;
           2.67415e-08 0.010425 0.00274946;
           2.13561e-09 0.00274946 0.0744749];
       
I3_world = T_03com(1:3,1:3) * I3_body * transpose(T_03com(1:3,1:3));
       
       
D = m1*transpose(Jvc1)*Jvc1 + m2*transpose(Jvc2)*Jvc2 + m3*transpose(Jvc3)*Jvc3  ...
    + transpose(Jw1)*I1_world*Jw1 + transpose(Jw2)*I2_world*Jw2 + transpose(Jw3)*I3_world*Jw3

syms theta1_dot theta2_dot theta3_dot theta4_dot pi
C_11 = calc_christoffel(D,1,1,1)*theta1_dot + calc_christoffel(D,2,1,1)*theta2_dot + calc_christoffel(D,3,1,1)*theta3_dot;
C_12 = calc_christoffel(D,1,2,1)*theta1_dot + calc_christoffel(D,2,2,1)*theta2_dot + calc_christoffel(D,3,2,1)*theta3_dot;
C_13 = calc_christoffel(D,1,3,1)*theta1_dot + calc_christoffel(D,2,3,1)*theta2_dot + calc_christoffel(D,3,3,1)*theta3_dot;

C_21 = calc_christoffel(D,1,1,2)*theta1_dot + calc_christoffel(D,2,1,2)*theta2_dot + calc_christoffel(D,3,1,2)*theta3_dot;
C_22 = calc_christoffel(D,1,2,2)*theta1_dot + calc_christoffel(D,2,2,2)*theta2_dot + calc_christoffel(D,3,2,2)*theta3_dot;
C_23 = calc_christoffel(D,1,3,2)*theta1_dot + calc_christoffel(D,2,3,2)*theta2_dot + calc_christoffel(D,3,3,2)*theta3_dot;

C_31 = calc_christoffel(D,1,1,3)*theta1_dot + calc_christoffel(D,2,1,3)*theta2_dot + calc_christoffel(D,3,1,3)*theta3_dot;
C_32 = calc_christoffel(D,1,2,3)*theta1_dot + calc_christoffel(D,2,2,3)*theta2_dot + calc_christoffel(D,3,2,3)*theta3_dot;
C_33 = calc_christoffel(D,1,3,3)*theta1_dot + calc_christoffel(D,2,3,3)*theta2_dot + calc_christoffel(D,3,3,3)*theta3_dot;

C = [C_11 C_12 C_13;
     C_21 C_22 C_23;
     C_31 C_32 C_33]

g = [0 0 9.81];

P1 = m1*g*T_01com(1:3,4) + m2*g*T_02com(1:3,4) + m3*g*T_03com(1:3,4);

G = [diff(P1,theta1);
     diff(P1,theta2);
     diff(P1,theta3);]
 
ht1 = matlabFunction(D)
ht2 = matlabFunction(C)
ht3 = matlabFunction(G)


%%

tau = zeros(3,1);

theta = zeros(3,1);
theta(2) = 0.1;
theta_dot = zeros(3,1);

theta_log = [];
theta_dot_log = [];

for i = 1:10000
    i
    tic
    [theta, theta_dot] = solve_dynamics(D,C,G,theta,theta_dot,tau);
    theta_log = [theta_log theta];
    theta_dot_log = [theta_dot_log theta_dot];
    toc
end



%%
function [C] = calc_christoffel(D,i,j,k)
    C = 0.5 * (diff(D(k,j),which_theta(i)) + diff(D(k,i),which_theta(j)) - diff(D(i,j),which_theta(k)));
end

function [theta] = which_theta(i)
    syms theta1 theta2 theta3 theta4
    if (i == 1)
        theta = theta1;
    elseif (i == 2)
        theta = theta2;
    elseif (i == 3)
        theta = theta3;
    end
end












