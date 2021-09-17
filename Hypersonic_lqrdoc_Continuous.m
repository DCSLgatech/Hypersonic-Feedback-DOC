function phaseout = Hypersonic_lqrdoc_Continuous(input)

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
rbar  = input.phase.state(:,1); % Radius
vbar  = input.phase.state(:,2); % Speed
fpa   = input.phase.state(:,3); % Flight path angle
dbar  = input.phase.state(:,4); % Downrange

tnum  = numel(rbar);

% u1dot = input.phase.control(:,1);

% Sensitivities
S1H   = input.phase.state(:,5);
S2H   = input.phase.state(:,6);
S3H   = input.phase.state(:,7);
S4H   = input.phase.state(:,8);

Svec = input.phase.state(:,5:8);
Pvec = input.phase.state(:, 9 : 18);

cbank = input.phase.control(:,1);

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

mass = input.auxdata.mass;

alpha    = input.auxdata.alpha;
SigmaP   = input.auxdata.SigmaP;
R        = input.auxdata.R;

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
df1du = zero;
df2du = zero;
df3du = -Lbar./(mass.*vbar);
df4du = zero;

% ---------------------------------------------------%
% --------- Evaluate Sensitivity Dynamics ---------- %
% ---------------------------------------------------%
SHdot = zeros(tnum,4);
PDOT  = zeros(tnum,10);

feedback_running_cost = zero; 



for i = 1:tnum
    
    % Get Riccati vector at time tk
    Pvec_k = Pvec(i, :);    
    % Reshape Riccati vector -> Riccati matrix
    P = [Pvec_k(1:4);
         0, Pvec_k(5:7);
         0, 0, Pvec_k(8:9);
         0, 0, 0, Pvec_k(10)];
    P = P + P'- diag([P(1,1), P(2,2), P(3,3), P(4,4)]);
%     P = reshape(Pvec_k, 4, 4);
    
    % Get sensitivity vector at time tk
    Svec_k = Svec(i, :);    
    % Reshape sensitivity vector -> sensitivity matrix
    S_k = reshape(Svec_k, 4, 1);
        
    S      = [S1H(i); S2H(i); S3H(i); S4H(i)];

    A      = [df1dx(i,:); df2dx(i,:); df3dx(i,:); df4dx(i,:)];
    B      = [df1dH(i); df2dH(i); df3dH(i); df4dH(i)];
    C      = [df1du(i); df2du(i); df3du(i); df4du(i)];
    K      = - R \ C' * P;
    Sdot   = (A + C*K)*S + B;
    Pdot   = -(A'*P + P*A - P*C/R*C'*P);
    
    SHdot(i,:) = Sdot';
%     PDOT(i,:)  = reshape(Pdot,1,16);
    PDOT(i,:)  = [Pdot(1,:), Pdot(2,2:4), Pdot(3,3:4), Pdot(4,4)];
    
    feedback_running_cost(i) = 0.5*trace(R * K * S * SigmaP * S' * K');
end

phaseout.dynamics  = [rbardot, vbardot, fpadot, dbardot, SHdot, PDOT];
phaseout.integrand = feedback_running_cost;
                  
end

    