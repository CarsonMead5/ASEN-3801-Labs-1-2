% Contributors: Victor Turpin Aguayo, Carson Mead, Liam Karlson, David Nguyen
% Course number: ASEN 3801
% File name: lab2
% Created: 1/20/26

clc;clear;close all

%% Question 2a - LoadASPENData and ConvertASPENData
% Reading in target and aerospace vehicle data and distributing into
% separate usable vectors

[t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData('3801_Sec001_test1.csv');

%% Question 2b - RotationMatrix321

% Testing
% attitude321 = [9 -2 33];
% DCM_321 = RotationMatrix321(attitude321)
% attitude321_2 = EulerAngles321(DCM_321)

%% Question 2c
% When beta (Y) is negative, it flips all the angles; when beta (Y) is 0,
% it just adds the Z angles (makes sense because no change)

% Testing
% attitude313 = [-9 0 -11]
% DCM_313 = RotationMatrix313(attitude313)
% attitude313_2 = EulerAngles313(DCM_313)

%% Question 5
% Orientation vs. Time

% Calculate 313 euler angles over time 
Euler_Angles313_ac = zeros(3,length(t_vec));
Euler_Angles313_tar = zeros(3, length(t_vec));
for i = 1:length(t_vec)
    Euler_Angles313_ac(:,i) = [EulerAngles313(RotationMatrix321(rad2deg(av_att(:,i))))];     % a/c 313 angles (alpha, beta, gamma)
    Euler_Angles313_tar(:,i) = [EulerAngles313(RotationMatrix321(rad2deg(tar_att(:,i))))];   % target 313 angles (alpha, beta, gamma)
end

%% Question 6-7
% Calculating Relative Position Vectors

% Calculating target relative to vehicle vector (Frame E)
pos_T_AV_inert = tar_pos_inert - av_pos_inert;

% Initializing relative position vector in body frame coordinates
pos_T_AV_body = zeros(size(pos_T_AV_inert));

% Converting position vector into Frame B coordinates
for i = 1:length(pos_T_AV_body)

    % Finding DCM to current attitude
    curr_DCM = RotationMatrix321(rad2deg(av_att(:,i)));

    % Finding current relative position vector
    pos_T_AV_body(:,i) = curr_DCM * pos_T_AV_inert(:,i);

end

%% Plots

% Question 3 (Frame N) Position vs. Time
figure();
plot3(av_pos_inert(1,:),av_pos_inert(2,:),av_pos_inert(3,:),'-b',LineWidth=1.25);
hold on;
plot3(tar_pos_inert(1,:),tar_pos_inert(2,:),tar_pos_inert(3,:),'--r',LineWidth=1.25);
grid on;
lgd = legend('Aerospace Vehicle','Target');
set(lgd,'Position',[0.6170157527116 0.304912513129444 0.271071432862963 0.0869047637212846]);
title("Vehicle and Target Paths")
xlabel('X Position (m)')
ylabel('Y Position (m)')
zlabel('Z Position (m)')
set(gca, 'ZDir', 'reverse');
print("3D Drone and Target Paths","-dpng","-r300")

% Question 4 (Frame E) Position and 3-2-1 Attitude vs. Time

% Position vs. Time
figure('Units','inches','Position',[2 2 9 4.5]);
sgtitle("Position vs. Time")
subplot(3, 1, 1);
hold on;
plot(t_vec, av_pos_inert(1,:),'-b',LineWidth=1.5);
plot(t_vec, tar_pos_inert(1,:),'--r',LineWidth=1.5);
hold off;
ylabel('X Position (m)')
lgd1 = legend('Aerospace Vehicle', 'Target');
set(lgd1,'Position',[0.661309520165125 0.911346132699368 0.24392857507297 0.0716666681198846]);
grid on;
ylim([-3,1.25])

subplot(3, 1, 2);
hold on;
plot(t_vec, av_pos_inert(2,:),'-b',LineWidth=1.5);
plot(t_vec, tar_pos_inert(2,:),'--r',LineWidth=1.5);
hold off;
ylabel('Y Position (m)')
grid on;
ylim([-7.5,0])

subplot(3, 1, 3);
hold on;
plot(t_vec, av_pos_inert(3,:),'-b',LineWidth=1.5);
plot(t_vec, tar_pos_inert(3,:),'--r',LineWidth=1.5);
hold off;
ylabel('Z Position (m)')
xlabel('Time(s)')
grid on;
ylim([-3,0.5])
saveas(gcf,"Position vs Time","png")

% 3-2-1 Euler Angles vs. Time
figure('Units','inches','Position',[2 2 9 4.5]);
sgtitle("3-2-1 Euler Angles vs. Time")
subplot(3, 1, 1);
hold on;
plot(t_vec, rad2deg(av_att(1,:)),'-b',LineWidth=1);
plot(t_vec, rad2deg(tar_att(1,:)),'--r',LineWidth=1);
hold off;
ylabel('Yaw \psi (deg)')
lgd2 = legend('Aerospace Vehicle', 'Target');
set(lgd2,'Position',[0.661309520165125 0.911346132699368 0.24392857507297 0.0716666681198846]);
grid on;

subplot(3, 1, 2);
hold on;
plot(t_vec, rad2deg(av_att(2,:)),'-b',LineWidth=1);
plot(t_vec, rad2deg(tar_att(2,:)),'--r',LineWidth=1);
hold off;
ylabel('Pitch \theta (deg)')
grid on;

subplot(3, 1, 3);
hold on;
plot(t_vec, rad2deg(av_att(3,:)),'-b',LineWidth=1);
plot(t_vec, rad2deg(tar_att(3,:)),'--r',LineWidth=1);
hold off;
ylabel('Roll \phi (deg)')
xlabel('Time(s)')
grid on;
saveas(gcf,'Euler 321 vs Time','png')

% Question 5 (3-1-3 Euler Angles over time)

figure('Units','inches','Position',[2 2 9 4.5]);
sgtitle("3-1-3 Euler Angles vs. Time")

subplot(3, 1, 1);
hold on;
plot(t_vec, Euler_Angles313_ac(1,:),'-b',LineWidth=1); % alpha (roll) of ac
plot(t_vec, Euler_Angles313_tar(1,:),'--r',LineWidth=1); % alpha (roll) of target
hold off;
ylabel('\alpha (deg)')
grid on;

subplot(3, 1, 2);
hold on;
plot(t_vec, Euler_Angles313_ac(2,:),'-b',LineWidth=1); % beta (pitch) of ac 
plot(t_vec, Euler_Angles313_tar(2,:),'--r',LineWidth=1); % beta (pitch) of target
hold off;
ylabel('\beta (deg)')
grid on;

subplot(3, 1, 3);
hold on;
plot(t_vec, Euler_Angles313_ac(3,:),'-b',LineWidth=1); % gamma (yaw) of ac
plot(t_vec, Euler_Angles313_tar(3,:),'--r',LineWidth=1); % gamma (yaw) of target
hold off;
ylabel('\gamma (deg)')
xlabel('Time (s)')
grid on;
lgd3 = legend('Aerospace Vehicle', 'Target');
set(lgd3,'Position',[0.724166663022268 0.90372708508032 0.24392857507297 0.0716666681198846]);
saveas(gcf,'Euler 313 vs Time','png')

% Question 6 (Frame E) Relative Position Vector vs. Time

figure('Units','inches','Position',[2 2 9 4.5])
sgtitle("Vehicle to Target Relative Position vs. Time")

subplot(3,1,1);
plot(t_vec,pos_T_AV_inert(1,:),'-b',LineWidth=1);
grid on;
ylabel("X Displacement (m)")
ylim([-3.5,3])

subplot(3,1,2);
plot(t_vec,pos_T_AV_inert(2,:),'-b',LineWidth=1);
grid on;
ylabel("Y Displacement (m)")
ylim([-5.5,3])

subplot(3,1,3);
plot(t_vec,pos_T_AV_inert(3,:),'-b',LineWidth=1);
grid on;
ylabel("Z Displacement (m)")
xlabel("Time (s)")
ylim([-2,1.25])
saveas(gcf,"Inertial Relative Position","png")

% Question 7 (Frame B) Relative Position Vector vs. Time

figure('Units','inches','Position',[2 2 9 4.5])
sgtitle('Target Relative Position in Vehicle Body Frame')

subplot(3,1,1)
plot(t_vec, pos_T_AV_body(1,:), '-b', LineWidth=1)
ylabel('X Displacement (m)'); 
grid on

subplot(3,1,2)
plot(t_vec, pos_T_AV_body(2,:), '-b', LineWidth=1.5)
ylabel('Y Displacement (m)'); 
grid on

subplot(3,1,3)
plot(t_vec, pos_T_AV_body(3,:), '-b', LineWidth=1.5)
ylabel('Z Displacement (m)'); 
xlabel('Time (s)'); 
grid on
saveas(gcf,"Body Relative Position","png")