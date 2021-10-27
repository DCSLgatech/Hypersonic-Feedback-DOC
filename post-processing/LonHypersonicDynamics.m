function xdot = LonHypersonicDynamics(t,x,ts,us,p,C)

% Extract state elements
rbar     = x(1);
phi      = x(2);
vbar     = x(3);
fpa      = x(4);

% Extract constants
rho0bar = C.rho0bar;
Sbar    = C.Sbar;
mubar   = C.mubar;
mass    = C.mass;
a0      = C.a0;
a1      = C.a1;
b0      = C.b0;
b1      = C.b1;
b2      = C.b2;

% Compute air density
rhobar = rho0bar * exp(-(rbar - 1) / p);

% Compute input at current time
alpha = interp1(ts, us, t);
alpha_hat = alpha * 180 / pi;

% Compute aerodynamic forces
CL   = a0 + a1 * alpha_hat;
CD   = b0 + b1 * alpha_hat + b2 * alpha_hat^2;
Dbar = 0.5 * rhobar * vbar^2 * Sbar * CD;
Lbar = 0.5 * rhobar * vbar^2 * Sbar * CL;

% Propagate dynamics
raddot = vbar .* sin(fpa);
phidot = vbar .* cos(fpa) ./ (rbar);
vdot   = -Dbar ./ mass - mubar .* sin(fpa) ./ rbar.^2;
fpadot = Lbar ./ (mass .* vbar) - mubar .* cos(fpa)./(rbar.^2 .* vbar) + vbar .* cos(fpa) ./ rbar;

% Output dynamics
xdot = [raddot; phidot; vdot; fpadot];

end