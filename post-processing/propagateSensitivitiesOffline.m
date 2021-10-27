function S = propagateSensitivitiesOffline(C, solution)

% Get state and control and time
state   = solution.phase.state;
rbar    = state(:, 1);
phi     = state(:, 2);
vbar    = state(:, 3);
fpa     = state(:, 4);

control   = solution.phase.control; 
alpha     = control;
alpha_hat = alpha * 180 / pi;

time      = solution.phase.time;
N       = length(time);

mubar     = C.mubar;
rho0bar   = C.rho0bar;
Hbar      = C.Hbar;
Sbar      = C.Sbar;

Re        = C.Re;
g0        = C.g0;

mass      = C.mass;
kQ        = C.kQ;

p         = C.Hbar;
R         = C.R;
SigmaP    = C.SigmaP;

a0        = C.a0;
a1        = C.a1;
b0        = C.b0;
b1        = C.b1;
b2        = C.b2;

% Some calculations
hbar      = (rbar - 1);
CD        = b0 + b1 * alpha_hat + b2 * alpha_hat.^2;
CL        = a0 + a1 * alpha_hat; 

rhobar    = rho0bar * exp(-hbar / Hbar);
qbar      = 0.5 * rhobar .* vbar .^ 2;
Dbar      = qbar .* Sbar .* CD;
Lbar      = qbar .* Sbar .* CL;

% Linearize dynamics
A_hist = zeros(N, 4, 4);
B_hist = zeros(N, 4, 1);
D_hist = zeros(N, 4, 1);

for k = 1 : N
    
    UQ_k = [Dbar(k), Lbar(k), mass, mubar, qbar(k), a1, b1, b2, hbar(k)]; % Useful quantities
    A_hist(k, :, :) = computeStateLinearization_v3(state(k,:)', control(k), p, UQ_k);
    B_hist(k, :, :) = computeControlLinearization_v3(state(k,:)', control(k), p, UQ_k);
    D_hist(k, :, :) = computeParameterLinearization_v3(state(k,:)', control(k), p, UQ_k);
    
end

% Propagate sensitivity dynamics
S0 = zeros(4, 1);
[~, S] = ode45(@(t, S) OLsensitivityDynamics(S, A_hist, D_hist, t, time), time, S0);

figure; 
subplot(1,4,1); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$\bar{S}_{r H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(1,4,2); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$\bar{S}_{\phi H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(1,4,3); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$\bar{S}_{v H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(1,4,4); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$\bar{S}_{\gamma H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(1,4,1);
plot(time, S(:, 1), '-k', 'LineWidth', 1);
%
subplot(1,4,2);
plot(time, S(:, 2), '-k', 'LineWIdth', 1);
%
subplot(1,4,3);
plot(time, S(:, 3), '-k', 'LineWidth', 1);
%
subplot(1,4,4);
plot(time, S(:, 4), '-k', 'LineWIdth', 1);


end