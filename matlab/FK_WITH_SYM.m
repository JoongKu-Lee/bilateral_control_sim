%% FK_WITH_SYM

clear all; clc;

syms theta1 theta2 theta3 
syms L01_z % L01_z = 0.07
syms L12_x L12_z % L12_x = 0.08, L12_z = 0.12
syms L23_y % L23_y = 0.43
syms L34_y L34_z % L34_y = 0.335 L34_z = -0.08
syms L1_xcom L1_zcom % L1_xcom = 0.008355, L1_zcom = 0.102659
syms L2_ycom L2_zcom % L2_ycom = 0.213502 L2_zcom = 0.0919281
syms L3_ycom L3_zcom % L3_ycom = 0.120975 L3_zcom = -0.0755858
syms m1 m2 m3 % m1 = 2.75458, m2 = 8.143, m3 = 5.14886
syms I1_xx I1_yy I1_zz I1_xy I1_yz I1_zx  % 0.0109697, 0.011856, 0.00604953, -2.08287e-08, 2.03974e-09, -0.000398309
syms I2_xx I2_yy I2_zz I2_xy I2_yz I2_zx  % 0.25374, 0.0212824, 0.247523, 1.41285e-07, 0.00499483, -7.13886e-08
syms I3_xx I3_yy I3_zz I3_xy I3_yz I3_zx  % 0.075434, 0.010425, 0.0744749, 2.67415e-08, 0.00274946, 2.13561e-09
syms g% 9.81

T_01 = [1 0 0 0;
        0 1 0 0;
        0 0 1 L01_z;
        0 0 0 1]; 

T_1dot2 = [0 0 1 L12_x;
           1 0 0 0;
           0 1 0 L12_z;
           0 0 0 1];

T_2dot3 = [1 0 0 0;
           0 1 0 L23_y;
           0 0 1 0;
           0 0 0 1];

T_3dot4 = [0 -1 0 0;
           0 0 1 L34_y;
           -1 0 0 L34_z;
           0 0 0 1];

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

T_01dot = T_01 * T_11dot;
T_02dot = T_01dot * T_1dot2 * T_22dot;
T_03dot = T_02dot * T_2dot3 * T_33dot;
T_0EE = T_03dot * T_3dot4;

T_1dot1com = [1 0 0 L1_xcom;
              0 1 0 0;
              0 0 1 L1_zcom;
              0 0 0 1];
T_01com = T_01dot * T_1dot1com;

T_2dot2com = [1 0 0 0;
              0 1 0 L2_ycom;
              0 0 1 L2_zcom;
              0 0 0 1];
T_02com = T_02dot * T_2dot2com;

T_3dot3com = [1 0 0 0;
              0 1 0 L3_ycom;
              0 0 1 L3_zcom;
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

I1_body = [I1_xx I1_xy I1_zx;
           I1_xy I1_yy I1_yz;
           I1_zx I1_yz I1_zz];
% I1_body = [I1_xx 0 0;
%            0 I1_yy 0;
%            0 0 I1_zz];

I1_world = T_01com(1:3,1:3) * I1_body * transpose(T_01com(1:3,1:3));

I2_body = [I2_xx I2_xy I2_zx;
           I2_xy I2_yy I2_yz;
           I2_zx I2_yz I2_zz];
% I2_body = [I2_xx 0 0;
%            0 I2_yy 0;
%            0 0 I2_zz];
       
I2_world = T_02com(1:3,1:3) * I2_body * transpose(T_02com(1:3,1:3));
       
I3_body = [I3_xx I3_xy I3_zx;
           I3_xy I3_yy I3_yz;
           I3_zx I3_yz I3_zz];
% I3_body = [I3_xx 0 0;
%            0 I3_yy 0;
%            0 0 I3_zz];
       
I3_world = T_03com(1:3,1:3) * I3_body * transpose(T_03com(1:3,1:3));
       
       
D = m1*transpose(Jvc1)*Jvc1 + m2*transpose(Jvc2)*Jvc2 + m3*transpose(Jvc3)*Jvc3  ...
    + transpose(Jw1)*I1_world*Jw1 + transpose(Jw2)*I2_world*Jw2 + transpose(Jw3)*I3_world*Jw3;

syms theta1_dot theta2_dot theta3_dot theta4_dot
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
     C_31 C_32 C_33];

gravity = [0 0 g];

P1 = m1*gravity*T_01com(1:3,4) + m2*gravity*T_02com(1:3,4) + m3*gravity*T_03com(1:3,4);

G = [diff(P1,theta1);
     diff(P1,theta2);
     diff(P1,theta3);];

ccode(D, 'file', 'Mass.txt'); 
ccode(C, 'file', 'Coriolis.txt');
ccode(G, 'file', 'Gravity.txt');

% D = simplify(D)
% C = simplify(C)
% G = simplify(G)

%% LINEAR PARAMETERIZE

pi1 = [m1 I1_xx I1_xy I1_zx I1_yy I1_yz I1_zz];
pi2 = [m2 I2_xx I2_xy I2_zx I2_yy I2_yz I2_zz];
pi3 = [m3 I3_xx I3_xy I3_zx I3_yy I3_yz I3_zz];
pi4 = [L01_z L12_x L12_z L23_y L34_y L34_z L1_xcom L1_zcom L2_ycom L2_zcom L3_ycom L3_zcom]

pi = transpose([pi1 pi2 pi3 pi4]);

syms theta1_ddot theta2_ddot theta3_ddot

targetEquation = D * [theta1_ddot; theta2_ddot; theta3_ddot] + C * [theta1_dot; theta2_dot; theta3_dot] + G;

[CC1, TT1] = coeffs(targetEquation(1), [pi])
[CC2, TT2] = coeffs(targetEquation(2), [pi])
[CC3, TT3] = coeffs(targetEquation(3), [pi])

A = union(TT1, TT2);
THETA = union(A, TT3);

for j = 1:length(THETA)
    idx1 = find(TT1 == THETA(j));
    idx2 = find(TT2 == THETA(j));
    idx3 = find(TT3 == THETA(j));
    
    if (length(idx1) == 0)
        Y(1,j) = 0;
    else
        Y(1,j) = CC1(idx1);
    end
    if (length(idx2) == 0)
        Y(2,j) = 0;
    else
        Y(2,j) = CC2(idx2);
    end
    if (length(idx3) == 0)
        Y(3,j) = 0;
    else
        Y(3,j) = CC3(idx3);
    end
end

%% LINEAR PARAMETERIZE FOR PASSIVITY BASED CONTROLLER

pi1 = [m1 I1_xx I1_xy I1_zx I1_yy I1_yz I1_zz];
pi2 = [m2 I2_xx I2_xy I2_zx I2_yy I2_yz I2_zz];
pi3 = [m3 I3_xx I3_xy I3_zx I3_yy I3_yz I3_zz];
pi4 = [L01_z L12_x L12_z L23_y L34_y L34_z L1_xcom L1_zcom L2_ycom L2_zcom L3_ycom L3_zcom];

pi = transpose([pi1 pi2 pi3 pi4]);

syms a1 a2 a3 v1 v2 v3

targetEquation = D * [a1; a2; a3] + C * [v1; v2; v3] + G;

[CC1, TT1] = coeffs(targetEquation(1), [pi])
[CC2, TT2] = coeffs(targetEquation(2), [pi])
[CC3, TT3] = coeffs(targetEquation(3), [pi])

A = union(TT1, TT2);
THETA = union(A, TT3);

for j = 1:length(THETA)
    idx1 = find(TT1 == THETA(j));
    idx2 = find(TT2 == THETA(j));
    idx3 = find(TT3 == THETA(j));
    
    if (length(idx1) == 0)
        Y(1,j) = 0;
    else
        Y(1,j) = CC1(idx1);
    end
    if (length(idx2) == 0)
        Y(2,j) = 0;
    else
        Y(2,j) = CC2(idx2);
    end
    if (length(idx3) == 0)
        Y(3,j) = 0;
    else
        Y(3,j) = CC3(idx3);
    end
end

Y 
THETA


%%

syms theta1 theta2 theta3 
syms L01_z % L01_z = 0.07
syms L12_x L12_z % L12_x = 0.08, L12_z = 0.12
syms L23_y % L23_y = 0.43
syms L34_y L34_z % L34_y = 0.335 L34_z = -0.08
syms L1_xcom L1_zcom % L1_xcom = 0.008355, L1_zcom = 0.102659
syms L2_ycom L2_zcom % L2_ycom = 0.213502 L2_zcom = 0.0919281
syms L3_ycom L3_zcom % L3_ycom = 0.120975 L3_zcom = -0.0755858
syms m1 m2 m3 % m1 = 2.75458, m2 = 8.143, m3 = 5.14886
syms I1_xx I1_yy I1_zz I1_xy I1_yz I1_zx  % 0.0109697, 0.011856, 0.00604953, -2.08287e-08, 2.03974e-09, -0.000398309
syms I2_xx I2_yy I2_zz I2_xy I2_yz I2_zx  % 0.25374, 0.0212824, 0.247523, 1.41285e-07, 0.00499483, -7.13886e-08
syms I3_xx I3_yy I3_zz I3_xy I3_yz I3_zx  % 0.075434, 0.010425, 0.0744749, 2.67415e-08, 0.00274946, 2.13561e-09

asd = subs(THETA, [L01_z, L12_x, L12_z, L23_y, L34_y, L34_z, L1_xcom, L1_zcom, L2_ycom, L2_zcom, L3_ycom, L3_zcom, m1, m2, m3, I1_xx, I1_yy, I1_zz, I1_xy, I1_yz, I1_zx, I2_xx, I2_yy, I2_zz, I2_xy, I2_yz, I2_zx, I3_xx, I3_yy, I3_zz, I3_xy, I3_yz, I3_zx], ...
    [0.07, 0.08, 0.12, 0.43, 0.335, -0.08, 0.008355, 0.102659, 0.213502, 0.0919281, 0.120975, -0.0755858, 2.75458, 8.143, 5.14886, 0.0109697, 0.011856, 0.00604953, -2.08287e-08, 2.03974e-09, -0.000398309, 0.25374, 0.0212824, 0.247523, 1.41285e-07, 0.00499483, -7.13886e-08, 0.075434, 0.010425, 0.0744749, 2.67415e-08, 0.00274946, 2.13561e-09])

double(asd)


%%

xp = T_0EE(1:3,4);
xp = simplify(xp);
Jv1 = diff(xp,theta1);
Jv2 = diff(xp,theta2);
Jv3 = diff(xp,theta3);
Jw1 = T_01dot(1:3,3);
Jw2 = T_02dot(1:3,3);
Jw3 = T_03dot(1:3,3);

Jacobian = [Jv1 Jv2 Jv3 ;
            Jw1 Jw2 Jw3 ];

Jacobian = simplify(Jacobian)

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






