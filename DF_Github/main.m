%%%%%%%%%%%%%%%%% plot curves
clear;clc;close all;
%% step1: Specify a smooth trajectory for flat ouputs and high derivatives
load('traj.mat');
subplot(411)
plot(traj.time,traj.x(1,:));grid on;title('ref\_xP[m]');
subplot(412)
plot(traj.time,traj.y(1,:));grid on;title('ref\_yP[m]');
subplot(413)
plot(traj.time,traj.z(1,:));grid on;title('ref\_zP[m]');
subplot(414)
plot(traj.time,180/pi*(traj.aB(1,:)));grid on;title('ref\_alpha2[degree]');
xlabel('time[s]');
%% Step2:  Generate X,U by using differential flatness
clc;close all;
[DF_out] = LG_DiffFlat(traj);% differential flatness function
AniDF = LG_AniDF(DF_out,8,0);% make a 3D animation