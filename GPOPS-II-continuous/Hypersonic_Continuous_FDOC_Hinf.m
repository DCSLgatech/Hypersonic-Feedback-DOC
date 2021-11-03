function phaseout = Hypersonic_Continuous_FDOC_Hinf(input)

% ---------------------------------------------------%
% ------ Extract Each Component of the State ------- %
% ---------------------------------------------------%
rbar        = input.phase.state(:,1);
phi         = input.phase.state(:,2);
vbar        = input.phase.state(:,3);
fpa         = input.phase.state(:,4);

gamma       = input.auxdata.gamma;

alpha       = input.phase.control;
alpha_hat   = alpha * 180 / pi;

Svec        = input.phase.state(:, 5 : 8);
Pvec        = input.phase.state(:, 9 : 24);

L           = length(rbar);

% ---------------------------------------------------%
% ----------------- Auxiliary data ------------------%
% ---------------------------------------------------%
mubar    = input.auxdata.mubar;
rho0bar  = input.auxdata.rho0bar;
Hbar     = input.auxdata.Hbar;
Sbar     = input.auxdata.Sbar;

Re       = input.auxdata.Re;
g0       = input.auxdata.g0;

mass     = input.auxdata.mass;
kQ       = input.auxdata.kQ;

p        = Hbar;
R        = input.auxdata.R;
Q        = input.auxdata.Q;
M        = input.auxdata.M;
SigmaP   = input.auxdata.SigmaP;

a0       = input.auxdata.a0;
a1       = input.auxdata.a1;
b0       = input.auxdata.b0;
b1       = input.auxdata.b1;
b2       = input.auxdata.b2;

% Some calculations
hbar     = (rbar - 1);
CD       = b0 + b1 * alpha_hat + b2 * alpha_hat.^2;
CL       = a0 + a1 * alpha_hat;

rhobar   = rho0bar * exp(-hbar / Hbar);
qbar     = 0.5 * rhobar .* vbar .^ 2;
Dbar     = qbar .* Sbar .* CD;
Lbar     = qbar .* Sbar .* CL;

% ---------------------------------------------------%
% ---- Evaluate Right-Hand Side of the Dynamics ---- %
% ---------------------------------------------------%

raddot = vbar .* sin(fpa);
phidot = vbar .* cos(fpa) ./ (rbar);
vdot   = -Dbar ./ mass - mubar .* sin(fpa) ./ rbar.^2;
fpadot = Lbar ./ (mass .* vbar) - mubar .* cos(fpa)./(rbar.^2 .* vbar) + vbar .* cos(fpa) ./ rbar;


% -----------------------------------------------------------------%
% ---------------------- Evaluate Feedback Cost ------------------ %
% -----------------------------------------------------------------%

% Initialize Riccati matrix dynamics
Pdot = zeros(L, 16);

% Initialize sensitivity matrix dynamics
Sdot = zeros(L, 4);

% Initialize feedback cost
J_feedback = zeros(L, 1);

for k = 1 : L % Loop through each time step 
    
    % Get Riccati vector at time tk
    Pvec_k = Pvec(k, :);
    
    % Reshape Riccati vector -> Riccati matrix
    P_k = reshape(Pvec_k, 4, 4);
    
    % Get sensitivity vector at time tk
    Svec_k = Svec(k, :);
    
    % Reshape sensitivity vector -> sensitivity matrix
    S_k = reshape(Svec_k, 4, 1);
    
    % Get state at time tk
    x_k = [rbar(k); phi(k); vbar(k); fpa(k)];
    
    % Get control at time tk
    u_k = alpha(k);
    
    % Compute linearizations at time tk
    UQ_k = [Dbar(k), Lbar(k), mass, mubar, qbar(k), a1, b1, b2, hbar(k), Sbar]; % Useful quantities
    A_k = computeStateLinearization_v3(x_k, u_k, p, UQ_k);
    B_k = computeControlLinearization_v3(x_k, u_k, p, UQ_k);
    D_k = computeParameterLinearization_v3(x_k, u_k, p, UQ_k);
    
    % Compute Kalman gain
    K_k = -R \ B_k' * P_k;
    L_k = (1 / gamma^2) * M^(-1) * D_k' * P_k;
    
    % Propagate sensitivity dynamics
    Sdot_k = (A_k + B_k * K_k) * S_k + D_k;
    
    % Reshape sensitivity matrix -> sensitivity vector
    Sdotvec_k = reshape(Sdot_k, 1, 4);
    
    % Propagate Riccati dynamics
    X_k    = B_k / R * B_k' - (1 / gamma^2) * (D_k * M^(-1) * D_k');
    Pdot_k = -(A_k' * P_k + P_k * A_k  - P_k * X_k * P_k + Q);

    % Reshape Riccati matrix -> Riccati vector
    Pdotvec_k = reshape(Pdot_k, 1, 16);
    
    % Store derivatives
    Sdot(k, :) = Sdotvec_k;
    Pdot(k, :) = Pdotvec_k;
    
    % Store running cost
    J_feedback(k) = trace((Q + K_k' * R * K_k - gamma^2 * L_k' * M * L_k) * S_k * SigmaP * S_k');

end

% ---------------------------------------------------%
% --------------- Path Constraints ----------------- %
% ---------------------------------------------------%

heating_rate     = kQ * sqrt(g0)^3 * sqrt(rhobar) .* vbar.^3 ;
dynamic_pressure = qbar .* g0 / Re^2;
normal_load      = sqrt(Lbar.^2 + Dbar.^2) * g0 / mass;

% Output dynamics, path constraints, and running cost
phaseout.dynamics  = [raddot, phidot, vdot, fpadot, Sdot, Pdot];
phaseout.path      = [heating_rate, dynamic_pressure, normal_load];
phaseout.integrand = J_feedback;

end

    