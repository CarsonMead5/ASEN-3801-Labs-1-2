% Author: Victor Turpin Aguayo
% Course number: ASEN 3801
% File name: RotationMatrix313
% Created: 1/24/26

function DCM = RotationMatrix313(attitude313)
%
% Inputs: attitude321 = 3x1 vector with the 3-1-3 Euler angles in the form 
%                       attitude313 = [alpha, beta, gamma]T in degrees
% 
% Outputs: DCM = 3x3 rotation matrix calculated from the Euler angles.
%
% Methodology: Use angle2dcm to obtain the DCM from the inputted Euler
% angles. As it is a 313 sequence, the order is 'ZXZ'.

attitude_rad = deg2rad(attitude313);
alpha = attitude_rad(1); beta = attitude_rad(2); gamma = attitude_rad(3);
DCM = angle2dcm(alpha, beta, gamma, 'ZXZ');

end