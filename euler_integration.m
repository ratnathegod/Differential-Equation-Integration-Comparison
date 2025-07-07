function [x, y] = euler_integration(g, init, xf, dx)
% EULER_INTEGRATION Numerically integrates an ODE using the explicit Euler method.
%
%   [x, y] = euler_integration(g, init, xf, dx)
%
%   Inputs:
%       g    - Function handle representing dy/dx = g(x,y)
%       init - Initial condition [x0, y0]
%       xf   - Final value of x
%       dx   - Step size for integration
%
%   Outputs:
%       x - Vector of x values from x0 to xf
%       y - Vector of approximated y values computed using Euler's method
%
%   Example:
%       g = @(x,y) -exp(-x.^2);
%       [x, y] = euler_integration(g, [0,1], 2, 0.1);

    x0 = init(1);
    y0 = init(2);
    x = x0:dx:xf;
    N = length(x);
    y = zeros(1, N);
    y(1) = y0;
    for i = 1:(N-1)
        y(i+1) = y(i) + g(x(i), y(i)) * dx;
    end
end
