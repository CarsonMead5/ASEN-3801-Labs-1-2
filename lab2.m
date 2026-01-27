% Contributors: Victor Turpin Aguayo, Carson Mead
% Course number: ASEN 3801
% File name: lab2
% Created: 1/20/26

clc;clear;close all

%% 2a - LoadASPENData and ConvertASPENData
% Reading in target and aerospace vehicle data and distributing into
% separate usable vectors

[t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData('3801_Sec001_test1.csv');

%% 2b - RotationMatrix321
% 

attitude321 = [9 -2 33];
DCM_321 = RotationMatrix321(attitude321)

attitude321_2 = EulerAngles321(DCM_321)
%% 2c
% When beta (Y) is negative, it flips all the angles; when beta (Y) is 0,
% it just adds the Z angles (makes sense because no change)
attitude313 = [-9 0 -11]
DCM_313 = RotationMatrix313(attitude313)

attitude313_2 = EulerAngles313(DCM_313)
%% 2d

%% 2e

%% Plots

% Question 3 (Frame N) Position vs. Time

M = readmatrix('3801_Sec001_test1.csv');
M(1:3,:) = []; % Delete header rows
M(:,2) = []; % Delete subframe column

pos_av_N = M(:,11:13)'/1000; % Convert to m (from mm)
pos_tar_N = M(:,5:7)'/1000; % Convert to m (from mm)

figure();
plot3(pos_av_N(1,:),pos_av_N(2,:),pos_av_N(3,:),'-b',LineWidth=1);
hold on
plot3(pos_tar_N(1,:),pos_tar_N(2,:),pos_tar_N(3,:),'--r',LineWidth=1);
grid on;
lgd = legend('Aerospace Vehicle','Target');
set(lgd,'Position',[0.374158609854457 0.756341084558016 0.271071432862963 0.0869047637212845]);
title("Vehicle and Target Paths")
xlabel('X position (m)')
ylabel('Y position (m)')
zlabel('Z position (m)')
print("3D Drone and Target Paths","-dpng","-r300")

% Question 4 (Frame E) Position and Attitude vs. Time

% Position vs. Time
figure();
sgtitle("Position vs. Time")
subplot(3, 1, 1);
plot(t_vec, av_pos_inert(1,:),'-b',LineWidth=1.5);
hold on;
plot(t_vec, tar_pos_inert(1,:),'--r',LineWidth=1.5);
ylabel('X position (m)')
lgd1 = legend('Aerospace Vehicle', 'Target');
set(lgd1,'Position',[0.661309520165125 0.911346132699368 0.24392857507297 0.0716666681198846]);
grid on;

subplot(3, 1, 2);
plot(t_vec, av_pos_inert(2,:),'-b',LineWidth=1.5);
hold on;
plot(t_vec, tar_pos_inert(2,:),'--r',LineWidth=1.5);
ylabel('Y position (m)')
grid on;

subplot(3, 1, 3);
plot(t_vec, av_pos_inert(3,:),'-b',LineWidth=1.5);
hold on;
plot(t_vec, tar_pos_inert(3,:),'--r',LineWidth=1.5);
ylabel('Z position (m)')
xlabel('Time(s)')
grid on;
saveas(gcf,"Position vs Time","png")

% Attitude vs. Time
figure();
sgtitle("Attitude vs. Time")
subplot(3, 1, 1);
plot(t_vec, rad2deg(av_att(3,:)),'-b',LineWidth=1);
hold on;
plot(t_vec, rad2deg(tar_att(3,:)),'--r',LineWidth=1);
ylabel('Yaw \psi (deg)')
lgd2 = legend('Aerospace Vehicle', 'Target');
set(lgd2,'Position',[0.661309520165125 0.911346132699368 0.24392857507297 0.0716666681198846]);
grid on;

subplot(3, 1, 2);
plot(t_vec, rad2deg(av_att(2,:)),'-b',LineWidth=1);
hold on;
plot(t_vec, rad2deg(tar_att(2,:)),'--r',LineWidth=1);
ylabel('Pitch \theta (deg)')
grid on;

subplot(3, 1, 3);
plot(t_vec, rad2deg(av_att(1,:)),'-b',LineWidth=1);
hold on;
plot(t_vec, rad2deg(tar_att(1,:)),'--r',LineWidth=1);
ylabel('Roll \phi (deg)')
xlabel('Time(s)')
grid on;
saveas(gcf,'Attitude vs Time','png')

% Question 5:
%calculate 313 euler angles over time 
Euler_Angles313_ac = zeros(3,length(t_vec));
Euler_Angles313_tar = zeros(3, length(t_vec));
for i = 1:length(t_vec)
    Euler_Angles313_ac(:,i) = [EulerAngles313(RotationMatrix321(rad2deg(av_att(:,i))))];     % a/c 313 angles (alpha, beta, gamma)
    Euler_Angles313_tar(:,i) = [EulerAngles313(RotationMatrix321(rad2deg(tar_att(:,i))))];   % target 313 angles (alpha, beta, gamma)
end

figure();
sgtitle("313 Euler Angles vs. Time")
subplot(3, 1, 1);
plot(t_vec, Euler_Angles313_ac(1,:),'-b',LineWidth=1); % alpha (roll) of ac
hold on; 
plot(t_vec, Euler_Angles313_tar(1,:),'--r',LineWidth=1); % alpha (roll) of target
ylabel('\alpha (deg)')
xlabel('Time(s)')
grid on;

subplot(3, 1, 2);
plot(t_vec, Euler_Angles313_ac(2,:),'-b',LineWidth=1); % beta (pitch) of ac
hold on; 
plot(t_vec, Euler_Angles313_tar(2,:),'--r',LineWidth=1); % beta (pitch) of target
ylabel('\beta (deg)')
xlabel('Time(s)')
grid on;

subplot(3, 1, 3);
title("313 Euler Angles for AC and Target")
plot(t_vec, Euler_Angles313_ac(3,:),'-b',LineWidth=1); % gamma (yaw) of ac
hold on; 
plot(t_vec, Euler_Angles313_tar(3,:),'--r',LineWidth=1); % gamma (yaw) of target
ylabel('\gamma (deg)')
xlabel('Time(s)')
grid on;
lgd3 = legend('Aerospace Vehicle', 'Target');
set(lgd3,'Position',[0.661309520165125 0.911346132699368 0.24392857507297 0.0716666681198846]);
saveas(gcf,'313EulerAngles','png')


% Question 6

