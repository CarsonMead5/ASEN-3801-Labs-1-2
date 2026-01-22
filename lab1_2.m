% Contributors: Victor Turpin Aguayo, Lola Feilke
% Course number: ASEN 3801
% File name: lab1_2.m
% Created: 1/20/26

clc;clear;close all

%% Establishing Key Problem Variables

% Initial state vector values
init.x = 0; % x-position (North)
init.y = 0; % y-position (East)
init.z = 0; % z-position (Down)
init.u = 0; % x-velocity (North)
init.v = 20; % y-velocity (East)
init.w = -20; % z-velocity (Down)

% Problem constants
const.Cd = 0.6; % Object coefficient of drag [n.d.]
const.diameter = 2; % Object diameter [cm]
const.A = pi*(const.diameter/100)^2; % Object cross-sectional area [m^2]
wind_vel = [0; 0; 0]; % Wind velocity [m/s]
const.m = 0.05; % Object mass [kg]
const.g = 9.81; % Gravitational acceleration [m/s^2]

% Calculate the density of the atmosphere in Boulder (@1655m)
[const.rho, const.a, const.temp, const.press] = stdatmo(1655);

% Time span
t_span = [0 20];

%% Running Runge-Kutta 4,5 Prediction with ode45

% Entering values into state vector
X_0 = [init.x; init.y; init.z; init.u; init.v; init.w];

% Adjusting tolerances and adding event function
options = odeset('Events',@termination_event,RelTol=1e-8,AbsTol=1e-8);

% Event function definition
function [value, isterminal, direction] = termination_event(~,X)
    value = (X(3) > 0); % Condition if it passes ground level (+z = down)
    isterminal = 1;   % Stop the integration
    direction  = 0;
end

% Running ode45
[t,X] = ode45(@(t,X) objectEOM(t,X,const.rho,const.Cd,const.A,const.m,const.g,wind_vel),t_span,X_0,options);

% Transposing resultant matrix (for reading convenience)
X = transpose(X);

%% Plotting Figures (Part C)

% Plotting Object Trajectory
figure()
plot3(X(1,:),X(2,:),X(3,:),LineWidth=1.5)
grid on
title("Object 3D Trajectory")
xlabel("Displacement (North)")
ylabel("Displacement (East)")
zlabel("Displacement (Down)")
set(gca, 'ZDir', 'reverse');
print(".\Lab1_2c_Figure","-dpng","-r300")

% Determining max height reached
max_height = max(abs(X(3,:)));

%% Varying Wind Vector (Part D)

% Vector of wind in the x-direction (North)
wind_x = 0:0.5:50;
wind_vel = zeros(3,length(wind_x));
wind_vel(1,:) = wind_x;

% Creating final position state vector matrix
final_pos = zeros(6,length(wind_x));

% Total final distance vector
final_distance = zeros(1,length(wind_x));

figure()
hold on

% Running ode45 for each wind condition
for i = 1:length(wind_x)

    % Running ode45
    [t_d,X_d] = ode45(@(t_d,X_d) objectEOM(t_d,X_d,const.rho,const.Cd,const.A,const.m,const.g,wind_vel(:,i)),t_span,X_0,options);
    X_d = transpose(X_d);

    % Storing final state
    final_pos(:,i) = X_d(:,end);
    
    % Storing final distance
    final_distance(i) = norm(final_pos(1:3,i));

    if (mod(i,10) == 0)
        plot(X_d(1,:),X_d(2,:),LineWidth=1.5);
    end

end

grid on
title("Object 3D Trajectory")
xlabel("Displacement (North)")
ylabel("Displacement (East)")
hold off

%% Plotting Figures (Part D)

% Plotting wind velocity effect on final x-position
figure()
plot(wind_x,final_pos(1,:),LineWidth=1.5)
grid on
title("Wind Velocity vs. Final X-Position")
ylabel("Horizontal Displacement (m)")
xlabel("Northbound Windspeed (m/s)")
print("./Crosswind vs X-Dist","-dpng","-r300")

% Plotting wind velocity effect on final distance from origin
figure()
plot(wind_x,final_distance,LineWidth=1.5)
grid on
title("Wind Velocity vs. Total Distance")
ylabel("Total Distance from Origin (m)")
xlabel("Northbound Windspeed (m/s)")
print("./Crosswind vs Total Dist","-dpng","-r300")