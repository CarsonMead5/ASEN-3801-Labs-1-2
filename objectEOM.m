% Contributors: Carson Mead
% Course Number: ASEN 3801
% File Name: objectEOM
% Created: 1/8/26

function xdot = objectEOM(t,x,rho,Cd,A,m,g,wind_vel)

% Inputs:
%   t - time (s)
%   x - state vector [position, velocity]
%       position - [x;y;z] (m) (object position (ICF))
%       velocity - [v_x;v_y;v_z] (object velocity (ICF)) (m/s)
%   rho - air density (kg/m^3)
%   Cd - coefficient of drag (n.d.)
%   A - cross-sectional area (m^2)
%   m - mass of object (kg)
%   g - gravitational acceleration (m/s^2)
%   wind_vel - wind velocity [w_x;w_y;w_z] (ICF)

% Outputs:
%   xdot - time derivative state vector [velocity,acceleration]
%       velocity - [v_x;v_y;v_z] (m/s)
%       acceleration - [a_z;a_y;a_z] (m/s^2)

% Extracting vector quantities from the state vector
position = x(1:3);
velocity = x(4:6);

% Calculating velocity of object relative to wind
rel_velocity = velocity - wind_vel;

% Calculating current airspeed (magnitude of velocity vector)
airspeed = norm(rel_velocity);

% Calculating magnitude of drag on object
D = 1/2*rho*airspeed^2*Cd*A;

% Calculating drag force vector
f_drag = -D / airspeed * rel_velocity;

% Establishing first-order ODES
v_x = velocity(1);
v_y = velocity(2);
v_z = velocity(3);
a_x = f_drag(1)/m;
a_y = f_drag(2)/m;
a_z = f_drag(3)/m + g;

% Building state vector xdot
xdot = [v_x;v_y;v_z;a_x;a_y;a_z];

end