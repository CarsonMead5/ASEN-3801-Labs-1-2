% Author: Victor Turpin Aguayo
% Course number: ASEN 3801
% File name: EulerAngles321
% Created: 1/24/26

function attitude321 = EulerAngles321(DCM)
%
% Inputs: DCM = 3x3 rotation matrix
% 
% Outputs: attitude321 = 3x1 vector with the 3-2-1 Euler angles in the form 
%                        attitude321 = [alpha, beta, gamma]T
%                        alpha = phi (roll angle)
%                        beta = theta (pitch angle)
%                        gamma = psi (yaw angle)
%
% Methodology: Use dcm2angle to obtain the Euler angles with the given
% rotation matrix. As it is a 321 sequence, we use the 'ZYX' rotation.
% However, we need to return the attitude in the [roll, pitch, yaw] order
% so we flip the array before returning it.

[yaw, pitch, roll] = dcm2angle(DCM,'ZYX');
attitude321 = rad2deg([roll; pitch; yaw]);

end