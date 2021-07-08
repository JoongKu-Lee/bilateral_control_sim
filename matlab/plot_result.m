%% PLOT RESULTS
close all
load('jointdata.mat')

%%
figure(1)
plot(out.simulink_joint1_log,'r:','LineWidth',2)
title('Joint1 Angle')
legend('By Simulation')

figure(2)
plot(out.simulink_joint2_log+0.1,'r:','LineWidth',2)
title('Joint2 Angle')
legend('By Simulation')

figure(3)
plot(out.simulink_joint3_log,'r:','LineWidth',2)
title('Joint3 Angle')
legend('By Simulation')

%%

figure(1)
plot([0.001:0.001:1],theta_log(1,1:1000),'b--','LineWidth',2)
title('Joint1 Angle')
legend('By Calculation')

figure(2)
plot([0.001:0.001:1],theta_log(2,1:1000),'b--','LineWidth',2)
title('Joint2 Angle')
legend('By Calculation')

figure(3)
plot([0.001:0.001:1],theta_log(3,1:1000),'b--','LineWidth',2)
title('Joint3 Angle')
legend('By Calculation')

%%
figure(1)
plot([0.001:0.001:1],theta_log(1,1:1000),'b--','LineWidth',2)
hold on
plot(out.simulink_joint1_log,'r:','LineWidth',2)
title('Joint1 Angle')
legend('By Calculation','By Simulation')

figure(2)
plot([0.001:0.001:1],theta_log(2,1:1000),'b--','LineWidth',2)
hold on
plot(out.simulink_joint2_log+0.1,'r:','LineWidth',2)
title('Joint2 Angle')
legend('By Calculation','By Simulation')

figure(3)
plot([0.001:0.001:1],theta_log(3,1:1000),'b--','LineWidth',2)
hold on
plot(out.simulink_joint3_log,'r:','LineWidth',2)
title('Joint3 Angle')
legend('By Calculation','By Simulation')


%%
clear all; clc;

a = load('jointdata.mat');
theta_log = a.theta_log;
theta_dot_log = a.theta_dot_log;

syms theta1 theta2 theta3

T_01 = SE3(SO3());
T_01.t = [0, 0, 0.07];

T_11dot = [cos(theta1) -sin(theta1) 0 0;
           sin(theta1) cos(theta1) 0 0;
           0 0 1 0
           0 0 0 1];

T_1dot2 = SE3(SO3([0, 0, 1; 1, 0, 0; 0, 1, 0]));
T_1dot2.t = [0.08, 0, 0.12];

T_22dot = [cos(theta2) -sin(theta2) 0 0;
           sin(theta2) cos(theta2) 0 0;
           0 0 1 0
           0 0 0 1];

T_2dot3 = SE3(SO3());
T_2dot3.t = [0, 0.43, 0];

T_33dot = [cos(theta3) -sin(theta3) 0 0;
           sin(theta3) cos(theta3) 0 0;
           0 0 1 0
           0 0 0 1];

T_3dot4 = SE3(SO3([1 0 0;0 1 0; 0 0 1]));
T_3dot4.t = [0, 0.335, -0.08];


T_01dot = T_01.T * T_11dot;
T_02 = T_01dot * T_1dot2.T;
T_02dot = T_02 * T_22dot;
T_03 = T_02dot * T_2dot3.T;
T_03dot = T_03 * T_33dot;
T_04 = T_03dot * T_3dot4.T;

figure(99)

for i = 1:10:4000
    joint1 = double(subs(T_01.T,[theta1, theta2, theta3],[theta_log(1,i), theta_log(2,i), theta_log(3,i)]));
    joint2 = double(subs(T_02,[theta1, theta2, theta3],[theta_log(1,i), theta_log(2,i), theta_log(3,i)]));
    joint3 = double(subs(T_03,[theta1, theta2, theta3],[theta_log(1,i), theta_log(2,i), theta_log(3,i)]));
    ee = double(subs(T_04,[theta1, theta2, theta3],[theta_log(1,i), theta_log(2,i), theta_log(3,i)]));

    joint1_pos = joint1(1:3,4)';
    joint2_pos = joint2(1:3,4)';
    joint3_pos = joint3(1:3,4)';
    ee_pos = ee(1:3,4)';
    
    X = [joint1_pos(1) joint2_pos(1) joint3_pos(1) ee_pos(1)]
    Y = [joint1_pos(2) joint2_pos(2) joint3_pos(2) ee_pos(2)]
    Z = [joint1_pos(3) joint2_pos(3) joint3_pos(3) ee_pos(3)]
    
    plot3(X,Y,Z,'-o','LineWidth',5,'Color','r')
    axis([-1 1 -1 1 -1 1])
    grid on
    pause(0.005)
end

