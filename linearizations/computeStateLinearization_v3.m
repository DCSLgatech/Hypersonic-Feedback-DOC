function A = computeStateLinearization_v3(x, u, p, UQ)

% State, control, and parameters
rbar  = x(1);
vbar  = x(3);
fpa   = x(4);
Hbar  = p(1);

% Useful quantities

Dbar   = UQ(1);
Lbar   = UQ(2);
mass   = UQ(3);
mubar  = UQ(4);

% Compute analytical Jacobian

% % % d (rdot) / d x % % %
df1dx = [0, 0, sin(fpa), vbar * cos(fpa)];

% % % d (phidot) / d x
df2dx = [-vbar * cos(fpa) / rbar^2, 0, cos(fpa) / rbar, -vbar * sin(fpa) / rbar];

% % % d (vdot) / d x
df3dx = [Dbar / (mass * Hbar) + 2 * mubar * sin(fpa) / (rbar ^ 3), 0, -2 * Dbar / (mass * vbar), 0];
     
% % % d (fpadot) / d x     
df4dx = [-Lbar  / (mass * vbar * Hbar) + 2 * mubar * cos(fpa) / (rbar ^3 * vbar) + ...
         -vbar * cos(fpa) / (rbar ^ 2), 0, Lbar / (mass * vbar ^ 2) + ...
         mubar * cos(fpa) / (rbar ^ 2 * vbar ^ 2) + cos(fpa) / rbar, ...
         mubar * sin(fpa) / (rbar ^ 2 * vbar) - vbar * sin(fpa) / rbar];
     
A     = [df1dx; df2dx; df3dx; df4dx];
     
end