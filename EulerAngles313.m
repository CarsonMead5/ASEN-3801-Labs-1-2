% Author: Victor Turpin Aguayo
% Course number: ASEN 3801
% File name: EulerAngles313
% Created: 1/24/26

function attitude313 = EulerAngles313(DCM)
%
% Inputs: DCM = 3x3 rotation matrix
% 
% Outputs: attitude321 = 3x1 vector with the 3-1-3 Euler angles in the form 
%                        attitude313 = [alpha, beta, gamma]T
%
% Methodology: Use dcm2angle to obtain the Euler angles with the given
% rotation matrix. As it is a 313 sequence, we use the 'ZXZ' rotation.

[alpha, beta, gamma] = dcm2angle(DCM,'ZXZ','Robust');
attitude313 = rad2deg([alpha; beta; gamma]);

end