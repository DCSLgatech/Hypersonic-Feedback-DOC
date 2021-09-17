function lqrsol = obtain_lqrgains(sol, data)

tspan = [sol.tau(end) sol.tau(1)];

x0 = reshape(data.Qf, 16, 1);

t_s     = sol.tau;
state_s = sol.xbar(:,1:4);
cbank_s = sol.cbank;

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[time, Pvec] = ode45(@(t,x)full_dynamics(t,x,t_s,state_s,cbank_s,data),...
                tspan,x0,options);
            
lqrsol.tau  = flip(time);
lqrsol.Pvec = flip(Pvec, 1) ;
lqrsol.K    = zeros(numel(time), 4);

for i = 1:numel(time)
    P = reshape(Pvec(i,:), 4, 4);
    eig(P);
    cond(P);
    zero = 0;
    
    % Sensitivities
    xbar = interp1(t_s, state_s, time(i));
    
    rbar = xbar(1);
    vbar = xbar(2);
    
    % ---------------------------------------------------%
    % ----------------- Auxiliary data ------------------%
    % ---------------------------------------------------%
    
    R       = data.R;
    rho0bar = data.rho0bar;
    Hbar    = data.Hbar;
    Sbar    = data.Sbar;
    
    CL   = data.CL;
    
    mass = data.mass;
    
    % Some calculations
    hbar = (rbar - 1);
    
    rhobar  = rho0bar * exp(-hbar / Hbar);
    
    qbar    = 0.5*rhobar.*vbar.^2;
    Lbar    = qbar.*Sbar.*CL;
    
    % bankdot = u1dot;
    
    % Partial of the dynamics with respect to control
    df1du = zero;
    df2du = zero;
    df3du = Lbar./(mass.*vbar);
    df4du = zero;
    
    C = [df1du; df2du; df3du; df4du];
    
    K = - inv(R) * C' * P;
    lqrsol.K(i,:) = K;
end

lqrsol.K = flip(lqrsol.K, 1);

end

function dxdt = full_dynamics(t, x, t_s, state_s, cbank_s, data)

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
P = reshape(x, 4, 4);

zero = 0;

% Sensitivities
xbar  = interp1(t_s, state_s, t);
cbank = interp1(t_s, cbank_s, t);

rbar  = xbar(1); % Radius
vbar  = xbar(2); % Speed
fpa   = xbar(3); % Flight path angle
dbar  = xbar(4); % Downrange

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

R    = data.R;

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
df1du = zero;
df2du = zero;
df3du = Lbar./(mass.*vbar);
df4du = zero;

% ---------------------------------------------------%
% --------- Evaluate Sensitivity Dynamics ---------- %
% ---------------------------------------------------%

A = [df1dx; df2dx; df3dx; df4dx];
C = [df1du; df2du; df3du; df4du];
Pdot = -(A'*P + P*A - P*C*inv(R)*C'*P);

dxdt = reshape(Pdot, 16, 1);


end