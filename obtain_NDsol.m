function NDsol = obtain_NDsol(solution, C, IC)

tspan = [solution.phase.time(1) solution.phase.time(end)]*ND.s2TU;

x0 = [IC.rbar; IC.vbar; IC.fpa; IC.dbar; zeros(4,1)];

t_s     = solution.phase.time;
cbank_s = solution.phase.control(:,1);

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[tbar, xbar] = ode45(@(t,x)full_dynamics(t,x,t_s,cbank_s,data),...
                tspan,x0,options);
            
sol.t    = tbar * ND.TU2s;         % Time
sol.SH   = [xbar(:,5), xbar(:,6)/ND.TU2s, xbar(:,7)/ND.TU2s, xbar(:,8)]; % Sensitivities

%% Process states for plotting
alt      = (xbar(:,1)*ND.DU2m - data.Rm)/1000;
speed    = xbar(:,2)*ND.DU2m/(ND.TU2s * 1000);
fpa      = xbar(:,3)*180/pi;
dran     = xbar(:,4)*ND.DU2m/1000;

sol.xplot = [alt, speed, fpa, dran];

end

%%
function dxdt = full_dynamics(t, x, t_s, cbank_s, data)

dxdt = zeros(8,1);

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
rbar  = x(1); % Radius
vbar  = x(2); % Speed
fpa   = x(3); % Flight path angle
dbar  = x(4); % Downrange

zero = 0;

% Sensitivities
S1H   = x(5);
S2H   = x(6);
S3H   = x(7);
S4H   = x(8);

cbank = interp1(t_s, cbank_s, t);

% ---------------------------------------------------%
% ----------------- Auxiliary data ------------------%
% ---------------------------------------------------%

mubar   = data.mubar;
rho0bar = data.rho0bar;
Hbar    = data.Hbar;
Sbar    = data.Sbar;

CD   = data.CD;
CL   = data.CL;
Rm   = data.Rm;
% P        = input.auxdata.P;
% Ra       = input.auxdata.Ra;
% Qf       = input.auxdata.Qf;
% Qa       = alpha * Qf;

mass = input.auxdata.mass;

% Some calculations
hbar = (rbar - 1);

rhobar  = rho0bar * exp(-hbar / Hbar);

qbar    = 0.5*rhobar.*vbar.^2;
Dbar    = qbar.*Sbar.*CD;
Lbar    = qbar.*Sbar.*CL;

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%
gbar = mubar./rbar.^2;

dxdt(1) = vbar.*sin(fpa);
dxdt(2) = -Dbar./mass - gbar.*sin(fpa);
dxdt(3) = Lbar.*cbank./(mass.*vbar) - gbar.*cos(fpa)./vbar + vbar.*cos(fpa)./rbar;
dxdt(4) = vbar.*cos(fpa);

% bankdot = u1dot;
%% Derivatives

% Partial of the dynamics with respect to state
df1dx = [zero, sin(fpa), vbar.*cos(fpa), zero];
df2dx = [(Dbar./(mass.*Hbar) + 2*gbar.*sin(fpa)./rbar), -2*Dbar./(mass.*vbar), ...
         -gbar.*cos(fpa), zero];
df3dx = [(-Lbar.*cbank./(mass.*vbar.*Hbar) + ...
         2*gbar.*cos(fpa)./(rbar.*vbar) - vbar.*cos(fpa)./(rbar.^2)),...
         (Lbar.*cbank./(mass.*vbar.^2) + gbar.*cos(fpa)./vbar.^2 + cos(fpa)./rbar),...
         (gbar.*sin(fpa)./vbar - vbar.*sin(fpa)./rbar), zero];
df4dx = [zero, cos(fpa), -vbar.*sin(fpa), zero];
     
% Partial of the dynamics with respect to uncertain parameters - H
df1dH = zero;
df2dH = -Dbar.*hbar./(mass.*Hbar.^2);
df3dH = Lbar.*cbank.*hbar./(mass.*vbar.*Hbar.^2);
df4dH = zero;

% Partial of the dynamics with respect to control
% df1du = zero;
% df2du = zero;
% df3du = -Lbar./(mass.*vbar);

% ---------------------------------------------------%
% --------- Evaluate Sensitivity Dynamics ---------- %
% ---------------------------------------------------%

S = [S1H; S2H; S3H; S4H];
A = [df1dx; df2dx; df3dx; df4dx];
B = [df1dH; df2dH; df3dH; df4dH];

Sdot = A * S + B;

dxdt(5:8,:) = Sdot;


end