function [theta, theta_dot] = solve_dynamics(D,C,G,theta,theta_dot,tau)
%SOLVE_DYNAMICS Calculates Dynamics

syms theta1 theta2 theta3 theta1_dot theta2_dot theta3_dot 

persistent theta_ddot;
if isempty(theta_ddot)
    theta_ddot = inv(D)*(tau-C*[theta1_dot; theta2_dot; theta3_dot]-G);
end

theta_ddot_calculated = subs(theta_ddot, [theta1, theta2, theta3, theta1_dot, theta2_dot, theta3_dot],[theta', theta_dot']);

theta_ddot_calculated = double(theta_ddot_calculated);

theta_dot = theta_dot + theta_ddot_calculated * 0.001;

theta = theta + theta_dot * 0.001;


end

%%