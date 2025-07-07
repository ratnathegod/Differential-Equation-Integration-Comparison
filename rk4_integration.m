function [x, y] = rk4_integration(g, init, xf, dx)
% RK4_INTEGRATION Numerically integrates an ODE using the 4th order Runge-Kutta method.
%
%   [x, y] = rk4_integration(g, init, xf, dx)
%
%   Inputs:
%       g    - Function handle representing dy/dx = g(x,y)
%       init - Initial condition [x0, y0]
%       xf   - Final value of x
%       dx   - Step size for integration
%
%   Outputs:
%       x - Vector of x values from x0 to xf
%       y - Vector of approximated y values computed using RK4
%
%   Example:
%       g = @(x,y) y*sin(x).^2;
%       [x, y] = rk4_integration(g, [0,pi], 3*pi, pi/8);

    x0 = init(1);
    y0 = init(2);
    x = x0:dx:xf;
    N = length(x);
    y = zeros(1, N);
    y(1) = y0;
    for i = 1:(N-1)
        k1 = g(x(i), y(i));
        k2 = g(x(i) + dx/2, y(i) + k1 * dx/2);
        k3 = g(x(i) + dx/2, y(i) + k2 * dx/2);
        k4 = g(x(i) + dx, y(i) + k3 * dx);
        y(i+1) = y(i) + (dx/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end
