function sol = reconstruct_solution(solution, data, ND, IC)

tspan = [solution.phase.time(1) solution.phase.time(end)];

x0 = [IC.rbar; IC.vbar; IC.fpa; IC.dbar; zeros(4,1)];

t_s     = solution.phase.time;
cbank_s = 1 * solution.phase.control(:,1);

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[time, state] = ode45(@(t,x)full_dynamics(t,x,t_s,cbank_s,data),...
                tspan,x0,options);
            
sol.tau     = time;         % Time
sol.xbar    = state(:,1:4); % States
sol.SHbar   = state(:,5:8); % Sensitivities

%% Process states for plotting
sol.tplot = sol.tau * ND.TU2s;

alt       = (state(:,1) * ND.DU2m - data.Rm)/1000;
speed     = state(:,2)* (ND.DU2m / ND.TU2s) / 1000;
fpa       = state(:,3)*180/pi;
dran      = state(:,4) * ND.DU2m / 1000;

sol.cbank = interp1(solution.phase.time,solution.phase.control(:,1),time);

sol.xplot  = [alt, speed, fpa, dran];
sol.SHplot = [state(:,5), state(:,6)/ND.TU2s, state(:,7)/ND.DU2m, state(:,8)]; 

end

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

mass = data.mass;

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

dxdt(1)  = vbar.*sin(fpa);
dxdt(2)  = -Dbar./mass - gbar.*sin(fpa);
dxdt(3)  = Lbar.*cbank./(mass.*vbar) - gbar.*cos(fpa)./vbar + vbar.*cos(fpa)./rbar;
dxdt(4)  = vbar.*cos(fpa);

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
%     C = [df1du(i); df2du(i); df3du(i)];
%     K = - Ra \ C' * M * S';
Sdot = A * S + B;

dxdt(5:8,:) = Sdot;


end