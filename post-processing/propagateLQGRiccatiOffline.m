function Pvec = propagateLQGRiccatiOffline(C, solution)

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
    
    UQ_k = [Dbar(k), Lbar(k), mass, mubar, qbar(k), a1, b1, b2, hbar(k), Sbar]; % Useful quantities
    A_hist(k, :, :) = computeStateLinearization_v3(state(k,:)', control(k), p, UQ_k);
    B_hist(k, :, :) = computeControlLinearization_v3(state(k,:)', control(k), p, UQ_k);
    D_hist(k, :, :) = computeParameterLinearization_v3(state(k,:)', control(k), p, UQ_k);
    
end

% Propagate sensitivity dynamics
Pf = C.Qf;
Pfvec = reshape(Pf, 16, 1);
[~, Pvec] = ode45(@(t, P) RiccatiDynamics(P, A_hist, B_hist, t, time, C), flipud(time), Pfvec);

% Plot time histories
figure; 
subplot(4,4,1); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{11}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,2); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{12}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,3); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{13}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,4); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{14}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,5); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{21}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,6); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{22}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,7); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{23}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,8); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{24}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,9); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{31}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,10); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{32}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,11); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{33}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,12); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{34}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,13); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{41}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,14); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{42}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,15); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{43}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,16); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$P_{44}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,4,1);
plot(flipud(time),Pvec(:,1),'-k','LineWidth',1);
%
subplot(4,4,2);
plot(flipud(time),Pvec(:,2),'-k','LineWidth',1);
%
subplot(4,4,3);
plot(flipud(time),Pvec(:,3),'-k','LineWidth',1);
%
subplot(4,4,4);
plot(flipud(time),Pvec(:,4),'-k','LineWidth',1);
%
subplot(4,4,5);
plot(flipud(time),Pvec(:,5),'-k','LineWidth',1);
%
subplot(4,4,6);
plot(flipud(time),Pvec(:,6),'-k','LineWidth',1);
%
subplot(4,4,7);
plot(flipud(time),Pvec(:,7),'-k','LineWidth',1);
%
subplot(4,4,8);
plot(flipud(time),Pvec(:,8),'-k','LineWidth',1);
%
subplot(4,4,9);
plot(flipud(time),Pvec(:,9),'-k','LineWidth',1);
%
subplot(4,4,10);
plot(flipud(time),Pvec(:,10),'-k','LineWidth',1);
%
subplot(4,4,11);
plot(flipud(time),Pvec(:,11),'-k','LineWidth',1);
%
subplot(4,4,12);
plot(flipud(time),Pvec(:,12),'-k','LineWidth',1);
%
subplot(4,4,13);
plot(flipud(time),Pvec(:,13),'-k','LineWidth',1);
%
subplot(4,4,14);
plot(flipud(time),Pvec(:,14),'-k','LineWidth',1);
%
subplot(4,4,15);
plot(flipud(time),Pvec(:,15),'-k','LineWidth',1);
%
subplot(4,4,16);
plot(flipud(time),Pvec(:,16),'-k','LineWidth',1);
%


end