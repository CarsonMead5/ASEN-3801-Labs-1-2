% Contributors: Carson Mead
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
%print(".\Lab1_2c_Figure","-dpng","-r300")

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

% Plotting Trajectories
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
    
    % Plotting current trajectory
    if (mod(wind_vel(1,i),10) == 0)
        plot3(X_d(1,:),X_d(2,:),X_d(3,:),LineWidth=1.5)
    end

end

hold off
grid on
title("Object Trajectories at Different Windspeeds")
xlabel("Displacement (North)")
ylabel("Displacement (East)")
zlabel("Displacement (Down)")
legend("0 m/s","10 m/s","20 m/s","30 m/s","40 m/s","50 m/s")
set(gca, 'ZDir', 'reverse');
view([-1.5,-1,0.5])
%print("./Object Varying Trajectories","-dpng","-r300")

%% Plotting Figures (Part D)

% Plotting wind velocity effect on final x-position
figure()
plot(wind_x,final_pos(1,:),LineWidth=1.5)
grid on
title("Wind Velocity vs. Final X-Position")
ylabel("Horizontal Displacement (m)")
xlabel("Northbound Windspeed (m/s)")
%print("./Crosswind vs X-Dist","-dpng","-r300")

% Plotting wind velocity effect on final distance from origin
figure()
plot(wind_x,final_distance,LineWidth=1.5)
grid on
title("Wind Velocity vs. Total Distance")
ylabel("Total Distance from Origin (m)")
xlabel("Northbound Windspeed (m/s)")
%print("./Crosswind vs Total Dist","-dpng","-r300")

%% Varying Geopotential Altitude and Wind Velocity (Part E)

% Creating array of altitudes and densities to test
altitudes = 0:100:4000;
rho = zeros(1,length(altitudes));

% Creating array to store minimum total distances
min_dist = zeros(1,length(altitudes));

% Populating different density values
for i = 1:length(rho)
    [rho(i),~,~,~] = stdatmo(altitudes(i));
end

% Running ode45 for each altitude and wind condition
figure() % Populating figure as ode45 is running
hold on

for i = 1:length(rho)
    for j = 1:length(wind_x)
    
        % Running ode45
        [t_e,X_e] = ode45(@(t_e,X_e) objectEOM(t_e,X_e,rho(i),const.Cd,const.A,const.m,const.g,wind_vel(:,j)),t_span,X_0,options);
        X_e = transpose(X_e);
    
        % Storing final state
        final_pos(:,j) = X_e(:,end);
        
        % Storing final distance
        final_distance(j) = norm(final_pos(1:3,j));
    
    end

    % Storing minimum total distance for a specific altitude
    [min_dist(i),min_dist_idx] = min(final_distance);

    % Plotting curve for current density
    if (mod(altitudes(i),800) == 0)
        h = plot(wind_x,final_distance,LineWidth=1.5);
        scatter(wind_vel(1,min_dist_idx),min_dist(i),20,h.Color,'filled')
    end
end
hold off
grid on
title("Wind Velocity vs. Total Distance at Different Altitudes")
xlabel("Northbound Windspeed (m/s)")
ylabel("Total Distance (m)")
legend("0 m","","800 m","","1600 m","","2400 m","","3200 m","","4000 m","","Location","northwest")
%print("./Total Distance at Different Altitudes","-dpng","-r300")

% Plotting minimum distance as a function of geopotential altitude
figure()
plot(altitudes,min_dist,LineWidth=1.5)
grid on
title("Minimum Distance vs. Geopotential Altitude")
xlabel("Geopotential Altitude (m)")
ylabel("Minimum Total Displacement (m)")
%print("./Minimum Total Distance at Different Altitudes","-dpng","-r300")

%% Varying Mass at a Constant Kinetic Energy (Part F)

% Determining kinetic energy value that must remain constant
KE = 1/2 * const.m * norm(X_0(4:6))^2;

% Varying mass
mass_vec = 0.005:0.0025:0.10; % [kg]

% Calculating initial velocity magnitude
V_0 = sqrt(2*KE./mass_vec);

% Creating a vector of initial states
X_0_vec = zeros(6,length(V_0));
X_0_vec(5,:) = sqrt(V_0.^2 / 2);
X_0_vec(6,:) = -sqrt(V_0.^2 / 2);

% Creating a new final distance vector
final_dist = zeros(1,length(V_0));

% Running ode45 for varied masses
for i = 1:length(mass_vec)

    % Running ode45
    [t_f,X_f] = ode45(@(t_f,X_f) objectEOM(t_f,X_f,const.rho,const.Cd,const.A,mass_vec(i),const.g,wind_vel(:,1)),t_span,X_0_vec(:,i),options);
    X_f = transpose(X_f);

    % Storing final state
    final_position = X_f(:,end);

    % Storing final distance
    final_dist(i) = norm(final_position);

end

% Plotting figure comparing final distance with mass
figure()
plot(mass_vec,final_dist,LineWidth=1.5)
grid on
title("Total Distance vs. Object Mass")
xlabel("Object Mass (kg)")
ylabel("Total Distance (m)")
%print("./Total Distance on Varied Masses","-dpng","-r300")

% Updating figure as ode45 runs
figure()
hold on

% Running ode45 for each initial velocity
for i = 1:length(V_0)
    for j = 1:length(wind_x) % Iterating through different wind conditions
    
        % Running ode45
        [t_f,X_f] = ode45(@(t_f,X_f) objectEOM(t_f,X_f,const.rho,const.Cd,const.A,mass_vec(i),const.g,wind_vel(:,j)),t_span,X_0_vec(:,i),options);
        X_f = transpose(X_f);
    
        % Storing final state
        final_pos(:,j) = X_f(:,end);
        
        % Storing final distance
        final_distance(j) = norm(final_pos(1:3,j));
    
    end

    % Plotting final distance vs wind velocity for specified mass
    if ((mod(mass_vec(i),0.005) == 0) && (mass_vec(i) < 0.045))
        plot(wind_vel(1,:),final_distance,LineWidth=1.5)
    end

end

hold off
grid on
title("Total Distance of Different Mass Objects at Varied Windspeeds")
xlabel("Northbound Windspeed (m/s)")
ylabel("Total Distance (m)")
legend("0.005 kg","0.01 kg","0.015 kg","0.02 kg","0.025 kg","0.03 kg","0.035 kg","0.04 kg","Location","northwest")
%print("./Varied Masses and Windspeeds","-dpng","-r300")