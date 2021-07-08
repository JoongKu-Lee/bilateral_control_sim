%% INITIALIZE ROBOT BY LOADING SIMULINK FILE
% BEFORE RUNNING, COMMENT OUT RIGIDBODYTREE RELATED THINGS
clear all; clc;

open_system('myrobot.slx')
[robot, importInfo] = importrobot(gcs)

