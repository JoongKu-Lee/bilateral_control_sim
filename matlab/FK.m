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

theta4 = 0;
T_44dot = SE3(SO3(rotz(theta4)));

T_01dot = T_01 * T_11dot;
T_02 = T_01dot * T_1dot2;
T_02dot = T_02 * T_22dot;
T_03 = T_02dot * T_2dot3;
T_03dot = T_03 * T_33dot;
T_04 = T_03dot * T_3dot4;
T_04dot = T_04 * T_44dot;

trplot(T_01, 'length', 0.2, 'color','k')
view(45, 45)
hold on
trplot(T_01dot, 'length', 0.2, 'color','b')
trplot(T_02 , 'length', 0.2, 'color','k')
trplot(T_02dot , 'length', 0.2, 'color','b')
trplot(T_03 , 'length', 0.2, 'color','k')
trplot(T_03dot , 'length', 0.2, 'color','b')
trplot(T_04 , 'length', 0.2, 'color','k')
trplot(T_04dot , 'length', 0.2, 'color','b')

%% JOINT JACOBIAN
clear all; clc;

syms theta1 theta2 theta3 theta4

T_01 = SE3(SO3());
T_01.t = [0, 0, 0.07];

T_1dot2 = SE3(SO3([0, 0, 1; 1, 0, 0; 0, 1, 0]));
T_1dot2.t = [0.08, 0, 0.12];

T_2dot3 = SE3(SO3());
T_2dot3.t = [0, 0.43, 0];

T_3dot4 = SE3(SO3([0, -1, 0; 0, 0, 1; -1, 0 ,0]));
T_3dot4.t = [0, 0.335, -0.08];

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
       
T_44dot = [cos(theta4) -sin(theta4) 0 0;
           sin(theta4) cos(theta4) 0 0;
           0 0 1 0
           0 0 0 1];

T_01dot = T_01.T * T_11dot;
T_02dot = T_01dot * T_1dot2.T * T_22dot;
T_03dot = T_02dot * T_2dot3.T * T_33dot;
T_04dot = T_03dot * T_3dot4.T * T_44dot;

xp = T_04dot(1:3,4);
Jv1 = diff(xp,theta1);
Jv2 = diff(xp,theta2);
Jv3 = diff(xp,theta3);
Jv4 = diff(xp,theta4);
Jw1 = T_01dot(1:3,3);
Jw2 = T_02dot(1:3,3);
Jw3 = T_03dot(1:3,3);
Jw4 = T_04dot(1:3,3);

Jacobian = [Jv1 Jv2 Jv3 Jv4;
            Jw1 Jw2 Jw3 Jw4]
        
%% COM JACOBIAN
clear all; clc;

syms theta1 theta2 theta3 theta4

T_01 = SE3(SO3());
T_01.t = [0, 0, 0.07];

T_1dot2 = SE3(SO3([0, 0, 1; 1, 0, 0; 0, 1, 0]));
T_1dot2.t = [0.08, 0, 0.12];

T_2dot3 = SE3(SO3());
T_2dot3.t = [0, 0.43, 0];

T_3dot4 = SE3(SO3([0, -1, 0; 0, 0, 1; -1, 0 ,0]));
T_3dot4.t = [0, 0.335, -0.08];

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
       
T_44dot = [cos(theta4) -sin(theta4) 0 0;
           sin(theta4) cos(theta4) 0 0;
           0 0 1 0
           0 0 0 1];

T_01dot = T_01.T * T_11dot;
T_02dot = T_01dot * T_1dot2.T * T_22dot;
T_03dot = T_02dot * T_2dot3.T * T_33dot;
T_04dot = T_03dot * T_3dot4.T * T_44dot;

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

T_4dot4com = [1 0 0 0;
              0 1 0 0;
              0 0 1 0.0316703;
              0 0 0 1];          
T_04com = T_04dot * T_4dot4com;
          
Pc1 = T_01com(1:3,4);
Jvc1 = [diff(Pc1, theta1) diff(Pc1, theta2) diff(Pc1, theta3) diff(Pc1, theta4)];
          
Pc2 = T_02com(1:3,4);
Jvc2 = [diff(Pc2, theta1) diff(Pc2, theta2) diff(Pc2, theta3) diff(Pc2, theta4)];

Pc3 = T_03com(1:3,4);
Jvc3 = [diff(Pc3, theta1) diff(Pc3, theta2) diff(Pc3, theta3) diff(Pc3, theta4)];

Pc4 = T_04com(1:3,4);
Jvc4 = [diff(Pc4, theta1) diff(Pc4, theta2) diff(Pc4, theta3) diff(Pc4, theta4)];
          
Jw1 = [T_01dot(1:3,3) zeros(3,1) zeros(3,1) zeros(3,1)];
Jw2 = [T_01dot(1:3,3) T_02dot(1:3,3) zeros(3,1) zeros(3,1)];
Jw3 = [T_01dot(1:3,3) T_02dot(1:3,3) T_03dot(1:3,3) zeros(3,1)];
Jw4 = [T_01dot(1:3,3) T_02dot(1:3,3) T_03dot(1:3,3) T_04dot(1:3,3)];

m1 = 2.75458;
m2 = 8.143;
m3 = 5.14886;
m4 = 0.749321;

I1_body = 1e-6 * [40009.404 0 -2760.658;
           0 41087.714 0;
           -2760.658 0 6244.520]; % 단위 kg mm^2

       
I1_world = T_01dot(1:3,1:3) * I1_body * transpose(T_01dot(1:3,1:3));

I2_body = 1e-6 * [694011.077 0 0;
           0 90142.949 -154892.125;
           0 -154892.125 618947.187]; % 단위 kg mm^2
       
I2_world = T_02dot(1:3,1:3) * I2_body * transpose(T_02dot(1:3,1:3));
       
I3_body = 1e-6 * [180280.4 0 0;
           0 39864.114 49849.954;
           0 49849.954 149891.072]; % 단위 kg mm^2
       
I3_world = T_03dot(1:3,1:3) * I3_body * transpose(T_03dot(1:3,1:3));
       
I4_body = 1e-6 * [1670.657 0 0;
           0 1670.657 0;
           0 0 1331.961]; % 단위 kg mm^2
       
I4_world = T_04dot(1:3,1:3) * I4_body * transpose(T_04dot(1:3,1:3));
       
D = m1*transpose(Jvc1)*Jvc1 + m2*transpose(Jvc2)*Jvc2 + m3*transpose(Jvc3)*Jvc3 + m4*transpose(Jvc4)*Jvc4 ...
    + transpose(Jw1)*I1_world*Jw1 + transpose(Jw2)*I2_world*Jw2 + transpose(Jw3)*I3_world*Jw3 + transpose(Jw4)*I4_world*Jw4;

syms theta1_dot theta2_dot theta3_dot theta4_dot
C_11 = calc_christoffel(D,1,1,1)*theta1_dot + calc_christoffel(D,2,1,1)*theta2_dot + calc_christoffel(D,3,1,1)*theta3_dot + calc_christoffel(D,4,1,1)*theta4_dot;
C_12 = calc_christoffel(D,1,2,1)*theta1_dot + calc_christoffel(D,2,2,1)*theta2_dot + calc_christoffel(D,3,2,1)*theta3_dot + calc_christoffel(D,4,2,1)*theta4_dot;
C_13 = calc_christoffel(D,1,3,1)*theta1_dot + calc_christoffel(D,2,3,1)*theta2_dot + calc_christoffel(D,3,3,1)*theta3_dot + calc_christoffel(D,4,3,1)*theta4_dot;
C_14 = calc_christoffel(D,1,4,1)*theta1_dot + calc_christoffel(D,2,4,1)*theta2_dot + calc_christoffel(D,3,4,1)*theta3_dot + calc_christoffel(D,4,4,1)*theta4_dot;

C_21 = calc_christoffel(D,1,1,2)*theta1_dot + calc_christoffel(D,2,1,2)*theta2_dot + calc_christoffel(D,3,1,2)*theta3_dot + calc_christoffel(D,4,1,2)*theta4_dot;
C_22 = calc_christoffel(D,1,2,2)*theta1_dot + calc_christoffel(D,2,2,2)*theta2_dot + calc_christoffel(D,3,2,2)*theta3_dot + calc_christoffel(D,4,2,2)*theta4_dot;
C_23 = calc_christoffel(D,1,3,2)*theta1_dot + calc_christoffel(D,2,3,2)*theta2_dot + calc_christoffel(D,3,3,2)*theta3_dot + calc_christoffel(D,4,3,2)*theta4_dot;
C_24 = calc_christoffel(D,1,4,2)*theta1_dot + calc_christoffel(D,2,4,2)*theta2_dot + calc_christoffel(D,3,4,2)*theta3_dot + calc_christoffel(D,4,4,2)*theta4_dot;

C_31 = calc_christoffel(D,1,1,3)*theta1_dot + calc_christoffel(D,2,1,3)*theta2_dot + calc_christoffel(D,3,1,3)*theta3_dot + calc_christoffel(D,4,1,3)*theta4_dot;
C_32 = calc_christoffel(D,1,2,3)*theta1_dot + calc_christoffel(D,2,2,3)*theta2_dot + calc_christoffel(D,3,2,3)*theta3_dot + calc_christoffel(D,4,2,3)*theta4_dot;
C_33 = calc_christoffel(D,1,3,3)*theta1_dot + calc_christoffel(D,2,3,3)*theta2_dot + calc_christoffel(D,3,3,3)*theta3_dot + calc_christoffel(D,4,3,3)*theta4_dot;
C_34 = calc_christoffel(D,1,4,3)*theta1_dot + calc_christoffel(D,2,4,3)*theta2_dot + calc_christoffel(D,3,4,3)*theta3_dot + calc_christoffel(D,4,4,3)*theta4_dot;

C_41 = calc_christoffel(D,1,1,4)*theta1_dot + calc_christoffel(D,2,1,4)*theta2_dot + calc_christoffel(D,3,1,4)*theta3_dot + calc_christoffel(D,4,1,4)*theta4_dot;
C_42 = calc_christoffel(D,1,2,4)*theta1_dot + calc_christoffel(D,2,2,4)*theta2_dot + calc_christoffel(D,3,2,4)*theta3_dot + calc_christoffel(D,4,2,4)*theta4_dot;
C_43 = calc_christoffel(D,1,3,4)*theta1_dot + calc_christoffel(D,2,3,4)*theta2_dot + calc_christoffel(D,3,3,4)*theta3_dot + calc_christoffel(D,4,3,4)*theta4_dot;
C_44 = calc_christoffel(D,1,4,4)*theta1_dot + calc_christoffel(D,2,4,4)*theta2_dot + calc_christoffel(D,3,4,4)*theta3_dot + calc_christoffel(D,4,4,4)*theta4_dot;

C = [C_11 C_12 C_13 C_14;
     C_21 C_22 C_23 C_24;
     C_31 C_32 C_33 C_34;
     C_41 C_42 C_43 C_44]

g = [0 0 -9.81];

P1 = m1*g*T_01com(1:3,4) + m2*g*T_02com(1:3,4) + m3*g*T_03com(1:3,4) + m4*g*T_04com(1:3,4);

G = [diff(P1,theta1);
     diff(P1,theta2);
     diff(P1,theta3);
     diff(P1,theta4)]



%%
function [C] = calc_christoffel(D,i,j,k)
    C = 0.5 * (diff(D(1,1),which_theta(i)) + diff(D(1,1),which_theta(j)) - diff(D(1,1),which_theta(k)));
end

function [theta] = which_theta(i)
    syms theta1 theta2 theta3 theta4
    if (i == 1)
        theta = theta1;
    elseif (i == 2)
        theta = theta2;
    elseif (i == 3)
        theta = theta3;
    elseif (i == 4)
        theta = theta4;
    end
end












