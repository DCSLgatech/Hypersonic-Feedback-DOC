function phaseout = Hypersonic_odoc_Continuous(input)

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
r     = input.phase.state(:,1); % Radius
v     = input.phase.state(:,2); % Speed
fpa   = input.phase.state(:,3); % Flight path angle
dran  = input.phase.state(:,4); % Downrange

tnum  = numel(r);

% u1dot = input.phase.control(:,1);

% Sensitivities
S1H   = input.phase.state(:,5);
S2H   = input.phase.state(:,6);
S3H   = input.phase.state(:,7);
S4H   = input.phase.state(:,8);

cbank = input.phase.control(:,1);

zero  = zeros(tnum,1);
unit  = ones(tnum,1);

% ---------------------------------------------------%
% ----------------- Auxiliary data ------------------%
% ---------------------------------------------------%

mu   = input.auxdata.mu;
rho0 = input.auxdata.rho0;
H    = input.auxdata.H;
S    = input.auxdata.S;

CD   = input.auxdata.CD;
CL   = input.auxdata.CL;
Rm   = input.auxdata.Rm;
% P        = input.auxdata.P;
% Ra       = input.auxdata.Ra;
% Qf       = input.auxdata.Qf;
% Qa       = alpha * Qf;

mass = input.auxdata.mass;

% Some calculations
h    = (r - Rm);

rho  = rho0 * exp(-h / H);

q    = 0.5*rho.*v.^2;
D    = q.*S.*CD;
L    = q.*S.*CL;

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%
g = mu./r.^2;

raddot  = v.*sin(fpa);
vdot    = -D./mass - g.*sin(fpa);
fpadot  = L.*cbank./(mass.*v) - g.*cos(fpa)./v + v.*cos(fpa)./r;
drandot = v.*cos(fpa);

% bankdot = u1dot;
%% Derivatives

% Partial of the dynamics with respect to state
df1dx = [zero, sin(fpa), v.*cos(fpa), zero];
df2dx = [(D./(mass.*H) + 2*g.*sin(fpa)./r), -2*D./(mass.*v), ...
         -g.*cos(fpa), zero];
df3dx = [(-L.*cbank./(mass.*v.*H) + ...
         2*g.*cos(fpa)./(r.*v) - v.*cos(fpa)./(r.^2)),...
         (L.*cbank./(mass.*v.^2) + g.*cos(fpa)./v.^2 + cos(fpa)./r),...
         (g.*sin(fpa)./v - v.*sin(fpa)./r), zero];
df4dx = [zero, cos(fpa), -v.*sin(fpa), zero];
     
% Partial of the dynamics with respect to uncertain parameters - H
df1dH = zero;
df2dH = -D.*h./(mass.*H.^2);
df3dH = L.*cbank.*h./(mass.*v.*H.^2);
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
    Sdot = A * S + B;
    
    SHdot(i,:) = Sdot';
end

phaseout.dynamics  = [raddot, vdot, fpadot, drandot, SHdot];
                  
end

    