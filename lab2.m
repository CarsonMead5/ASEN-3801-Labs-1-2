% Contributors: Victor Turpin Aguayo
% Course number: ASEN 3801
% File name: lab2
% Created: 1/20/26

clc;clear;close all

%% 2a
[t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData('3801_Sec001_test1.csv');

%% 2b

%% 2c

%% 2d

%% 2e

%% Plots
% Question 3 (Frame N)

M = readmatrix('3801_Sec001_test1.csv');
M(1:3,:) = []; % Delete header rows
M(:,2) = []; % Delete subframe column

pos_av_N = M(:,11:13)'/1000; % Convert to m (from mm)
pos_tar_N = M(:,5:7)'/1000; % Convert to m (from mm)

figure();
plot3(pos_av_N(1,:),pos_av_N(2,:),pos_av_N(3,:),'-b');
hold on
plot3(pos_tar_N(1,:),pos_tar_N(2,:),pos_tar_N(3,:),'--r');
grid on;
legend('Aerospace Vehicle','Target')
xlabel('X position (m)')
ylabel('Y position (m)')
zlabel('Z position (m)')

% Question 4 (Frame E)

figure();
subplot(3, 1, 1);
plot(t_vec, av_pos_inert(1,:),'-b');
hold on;
plot(t_vec, tar_pos_inert(1,:),'--r');
ylabel('X position (m)')
legend('Aerospace Vehicle', 'Target',Location='bestoutside')
grid on;

subplot(3, 1, 2);
plot(t_vec, av_pos_inert(2,:),'-b');
hold on;
plot(t_vec, tar_pos_inert(2,:),'--r');
ylabel('Y position (m)')
grid on;

subplot(3, 1, 3);
plot(t_vec, av_pos_inert(3,:),'-b');
hold on;
plot(t_vec, tar_pos_inert(3,:),'--r');
ylabel('Z position (m)')
xlabel('Time(s)')
grid on;

%% Need to change to degrees
figure();
subplot(3, 1, 1);
plot(t_vec, rad2deg(av_att(1,:)),'-b');
hold on;
plot(t_vec, rad2deg(tar_att(1,:)),'--r');
ylabel('Yaw \psi (deg)')
legend('Aerospace Vehicle', 'Target',Location='bestoutside')
grid on;

subplot(3, 1, 2);
plot(t_vec, rad2deg(av_att(2,:)),'-b');
hold on;
plot(t_vec, rad2deg(tar_att(2,:)),'--r');
ylabel('Pitch \theta (deg)')
grid on;

subplot(3, 1, 3);
plot(t_vec, rad2deg(av_att(3,:)),'-b');
hold on;
plot(t_vec, rad2deg(tar_att(3,:)),'--r');
ylabel('Roll \phi (deg)')
xlabel('Time(s)')
grid on;
