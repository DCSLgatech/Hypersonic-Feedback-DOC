function phaseout = Hypersonic_baseline_Continuous(input)

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
rbar  = input.phase.state(:,1); % Radius
vbar  = input.phase.state(:,2); % Speed
fpa   = input.phase.state(:,3); % Flight path angle
dbar  = input.phase.state(:,4); % Downrange

% Sensitivities
S1H   = input.phase.state(:,5);
S2H   = input.phase.state(:,6);
S3H   = input.phase.state(:,7);
S4H   = input.phase.state(:,8);

% u1dot = input.phase.control(:,1);

cbank = input.phase.control(:,1);

tnum  = numel(rbar);

zero  = zeros(tnum,1);
unit  = ones(tnum,1);

% ---------------------------------------------------%
% ----------------- Auxiliary data ------------------%
% ---------------------------------------------------%

mubar   = input.auxdata.mubar;
rho0bar = input.auxdata.rho0bar;
Hbar    = input.auxdata.Hbar;
Sbar    = input.auxdata.Sbar;

CD   = input.auxdata.CD;
CL   = input.auxdata.CL;
Rm   = input.auxdata.Rm;
% P        = input.auxdata.P;
% Ra       = input.auxdata.Ra;
% Qf       = input.auxdata.Qf;
% Qa       = alpha * Qf;

mass     = input.auxdata.mass;

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

rbardot  = vbar.*sin(fpa);
vbardot  = -Dbar./mass - gbar.*sin(fpa);
fpadot   = Lbar.*cbank./(mass.*vbar) - gbar.*cos(fpa)./vbar + vbar.*cos(fpa)./rbar;
dbardot  = vbar.*cos(fpa);

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
SHdot = zeros(tnum,4);

% running_cost          = zero;
% feedback_running_cost = zero;

for i = 1:tnum
    S = [S1H(i); S2H(i); S3H(i); S4H(i)];
    A = [df1dx(i,:); df2dx(i,:); df3dx(i,:); df4dx(i,:)];
    B = [df1dH(i); df2dH(i); df3dH(i); df4dH(i)];
%     C = [df1du(i); df2du(i); df3du(i)];
%     K = - Ra \ C' * M * S';
    Sdot = A * S + B;
    
    SHdot(i,:) = Sdot';
    
%     running_cost(i) = 0;
%     feedback_running_cost(i) = 0.5*trace(K'*kron(eye(1),Ra)*K);
end

phaseout.dynamics  = [rbardot, vbardot, fpadot, dbardot, SHdot];
% phaseout.integrand = running_cost + feedback_running_cost;
                  
end

    