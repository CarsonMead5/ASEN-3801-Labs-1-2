% Author: Victor Turpin Aguayo
% Course number: ASEN 3801
% File name: RotationMatrix321
% Created: 1/20/26

function DCM = RotationMatrix321(attitude321)
%
% Inputs: attitude321 = 3x1 vector with the 3-2-1 Euler angles in the form 
%                       attitude321 = [alpha, beta, gamma]T in degrees
%                       alpha = phi (roll angle)
%                       beta = theta (pitch angle)
%                       gamma = psi (yaw angle)
% 
% Outputs: DCM = 3x3 rotation matrix calculated from the Euler angles.
%
% Methodology: Use angle2dcm to obtain the DCM from the inputted Euler
% angles. As the attitude is expressed in the [roll, pitch, yaw], and not 
% the common [yaw, pitch, roll], we flip the array to have yaw first. The 
% 321 sequence has to be 'ZYX'.

attitude_rad = deg2rad(attitude321);
roll = attitude_rad(1); pitch = attitude_rad(2); yaw = attitude_rad(3); 
DCM = angle2dcm(yaw, pitch, roll, 'ZYX');

end