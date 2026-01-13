function d_statevector_dt = odefun(~,statevector)
    % Inputs: statevector = column vector of the state variables
    % = [w,x,y,z]'
    %
    % Outputs: d_statevector_dt = [w_dot,x_dot,y_dot,z_dot]'
    %
    % Methodology: Use the given differential equations to calculate the
    % statevector variables

    w = statevector(1); x = statevector(2); y = statevector(3); z = statevector(4);

    w_dot = -9*w + y;
    x_dot = 4*w*x*y - x^2;
    y_dot = 2*w - x - 2*z;
    z_dot = x*y - y^2 - 3*z^3;
    d_statevector_dt = [w_dot; x_dot; y_dot; z_dot];
end