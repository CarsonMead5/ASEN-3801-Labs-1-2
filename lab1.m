% Contributors: Victor Turpin Aguayo, Lola Feilke
% Course number: ASEN 3801
% File name: lab1
% Created: 1/12/26

clc;clear;close all

%% 1a
% Set initial conditions and time span
w0 = 1; x0 = 3; y0 = 2; z0 = 6;
statevector_0 = [w0 x0 y0 z0]';
tspan = [0 20];

options = odeset(RelTol=1e-8,AbsTol=1e-8);
[t,statevector] = ode45(@odefun,tspan,statevector_0,options);

w = statevector(:,1);
x = statevector(:,2);
y = statevector(:,3);
z = statevector(:,4);

figure(); 

% 1. Plot w (Top)
subplot(4, 1, 1);
plot(t, w, 'LineWidth', 1.5, 'Color', '#0072BD'); % Blue
ylabel('w (n.d.)');
title('Problem 1a: Time evolution of w, x, y, and z'); % Title only on the top plot
grid on;
set(gca, 'XTickLabel', []); % Hide x-axis labels for neat stacking
xlim([min(t) max(t)]);      % Tighten x-axis

% 2. Plot x
subplot(4, 1, 2);
plot(t, x, 'LineWidth', 1.5, 'Color', '#D95319'); % Orange
ylabel('x (n.d.)');
grid on;
set(gca, 'XTickLabel', []);
xlim([min(t) max(t)]);

% 3. Plot y
subplot(4, 1, 3);
plot(t, y, 'LineWidth', 1.5, 'Color', '#EDB120'); % Yellow/Gold
ylabel('y (n.d.)');
grid on;
set(gca, 'XTickLabel', []);
xlim([min(t) max(t)]);

% 4. Plot z (Bottom)
subplot(4, 1, 4);
plot(t, z, 'LineWidth', 1.5, 'Color', '#7E2F8E'); % Purple
ylabel('z (n.d.)');
xlabel('t (n.d.)'); % Label x-axis only at the bottom
grid on;
xlim([min(t) max(t)]);

%% 1b
% Calculate the reference values
tol_r = 1e-12;
options = odeset(RelTol=tol_r,AbsTol=tol_r);
[t,statevector] = ode45(@odefun,tspan,statevector_0,options);

w_r = statevector(end,1);
x_r = statevector(end,2);
y_r = statevector(end,3);
z_r = statevector(end,4);

tol = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10];
diff = zeros(4,length(tol));

% Get the difference for each different tolerance
for i = 1:length(tol)
options = odeset(RelTol=tol(i),AbsTol=tol(i));
[t,statevector] = ode45(@odefun,tspan,statevector_0,options);

w_end = statevector(end,1);
x_end = statevector(end,2);
y_end = statevector(end,3);
z_end = statevector(end,4);

diff(:,i) = abs([w_end - w_r; x_end - x_r; y_end - y_r; z_end - z_r]);
end

% Create the table (each row is a statevector variable)
variableNames = {'1e-2', '1e-4', '1e-6', '1e-8', '1e-10'};
T = array2table(diff, 'VariableNames', variableNames);

disp(T)

%% 2b
% Calculate the density of the atmosphere in Boulder (@1655m)
[rho, a, temp, press] = stdatmo(1655);
