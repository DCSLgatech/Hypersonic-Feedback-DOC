function closedloopS = obtain_closedloopS(data, ND, sol, feedback)

tspan = [sol.tau(1) sol.tau(end)];

x0 = zeros(4,1);

t_s     = sol.tau;
state_s = sol.xbar;
cbank_s = sol.cbank;

t_K     = feedback.tau;
K_s     = feedback.K;

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[time, state] = ode45(@(t,x)full_dynamics(t,x,t_s,state_s, cbank_s, t_K,...
        K_s, data),tspan,x0,options);
    
closedloopS.t  = time * ND.TU2s;
closedloopS.SH = [state(:,1), state(:,2)/ND.TU2s, state(:,3)/ND.DU2m, state(:,4)]; 

K = interp1(t_K, K_s, time);

KS = sum(K.*closedloopS.SH, 2);

closedloopS.CC = data.SigmaP * KS.^2;

end


%%
function dxdt = full_dynamics(t,x,t_s,state_s, cbank_s, t_K, K_s, data)

xbar = interp1(t_s, state_s, t);

K   = interp1(t_K, K_s, t);

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
rbar  = xbar(1); % Radius
vbar  = xbar(2); % Speed
fpa   = xbar(3); % Flight path angle
dbar  = xbar(4); % Downrange

cbank = interp1(t_s, cbank_s, t);

zero = 0;

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
     
% Partial of the dynamics with respect to uncertain parameters - H
df1du = zero;
df2du = zero;
df3du = Lbar./(mass.*vbar);
df4du = zero;

% ---------------------------------------------------%
% --------- Evaluate Sensitivity Dynamics ---------- %
% ---------------------------------------------------%

SH = [x(1); x(2); x(3); x(4)];
A = [df1dx; df2dx; df3dx; df4dx];
B = [df1dH; df2dH; df3dH; df4dH];
C = [df1du; df2du; df3du; df4du];
dxdt = (A + C*K)*SH + B;

end