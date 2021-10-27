function [muf, Sigmaf] = MonteCarlo(solution,C,ND,IC,col)

%% Monte Carlo Analysis
%close all;

% Load in solutions
%OL_solution = load('solution_OL_a1.mat');

% Generate 100 samples for parameters
num_MC = 100;
p_MC = normrnd(C.Hbar,sqrt(C.SigmaP),num_MC,1);

% Plot samples
figure(1); hold on; grid on; box on;
xlabel('Scale Height (km)','FontSize',15);
set(gca,'FontSize',15);
scatter(p_MC * ND.DU2m / 1000, zeros(num_MC,1), 'k+');

% Get open loop optimal control
tvec  = solution.phase.time;
u     = solution.phase.control;

% Construct time vector
M = 100;
tspan = linspace(0, tvec(end), M);

% Initial state
rbar0  = IC.rbar;
phi0   = IC.phi;
vbar0  = IC.vbar;
gamma0 = IC.fpa;
x0bar  = [rbar0; phi0; vbar0; gamma0];

% Initialize MC state histories
xhist = zeros(M, 4, num_MC);

% ODE45 options
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Loop through each parameter value
for ii = 1 : num_MC
    
    % Get MC parameter value
    p_ii = p_MC(ii);
    
    % Propagate dynamics 
    [~, xout]   = ode45(@(t,x) LonHypersonicDynamics(t,x,tvec,u,p_ii,C), tspan, x0bar, options);

    % Store MC trajectory
    xhist(:, :, ii) = xout;
    
end

% Interpolate nominal trajectories
x         = solution.phase.state;
x_nominal = zeros(M, 4);
u_nominal = zeros(M, 1);
%Svec         = zeros(M, 4);
for k = 1 : M
    
    tk = tspan(k);
    
    full_state_k    = interp1(tvec, x, tk);
    x_nominal(k, :) = full_state_k(1:4);
    u_nominal(k)    = interp1(tvec, u, tk);
    %Svec(k ,:)      = full_state_k(5:8);
    
end

% figure(2); 
% subplot(2,2,1); hold on; grid on;
% plot(tspan*ND.TU2s/60,Svec(:,1),'-k','LineWidth',1);
% subplot(2,2,2); hold on; grid on;
% plot(tspan*ND.TU2s/60,Svec(:,2),'-k','LineWidth',1);
% subplot(2,2,3); hold on; grid on;
% plot(tspan*ND.TU2s/60,Svec(:,3),'-k','LineWidth',1);
% subplot(2,2,4); hold on; grid on;
% plot(tspan*ND.TU2s/60,Svec(:,4),'-k','LineWidth',1);

% Plot nominal and MC trajectories
grey = [0.7 0.7 0.7];

figure(3); 
subplot(2,2,1); hold on; grid on;
xlabel('Time (min)','FontSize',15);
ylabel('Height (km)','FontSize',15);
set(gca,'FontSize',15);
subplot(2,2,2); hold on; grid on;
xlabel('Time (min)','FontSize',15);
ylabel('Longitude (deg)','FontSize',15);
set(gca,'FontSize',15);
subplot(2,2,3); hold on; grid on;
xlabel('Time (min)','FontSize',15);
ylabel('Speed (km/s)','FontSize',15);
set(gca,'FontSize',15);
subplot(2,2,4); hold on; grid on;
xlabel('Time (min)','FontSize',15);
ylabel('Flight path angle (deg)','FontSize',15);
set(gca,'FontSize',15);

subplot(2,2,1);
for ii = 1 : num_MC
    plot(tspan*ND.TU2s/60,(xhist(:,1,ii)-1)*ND.DU2m/1000,'color',grey);
end
plot(tspan*ND.TU2s/60,(x_nominal(:,1)-1)*ND.DU2m/1000,'color',col,'LineWidth',1);
%
subplot(2,2,2);
for ii = 1 : num_MC
    plot(tspan*ND.TU2s/60,xhist(:,2,ii)*180/pi,'color',grey);
end
plot(tspan*ND.TU2s/60,x_nominal(:,2)*180/pi,'color',col,'LineWidth',1);
%
subplot(2,2,3);
for ii = 1 : num_MC
    plot(tspan*ND.TU2s/60,xhist(:,3,ii)*ND.DU2m/ND.TU2s/1000,'color',grey);
end
plot(tspan*ND.TU2s/60,x_nominal(:,3)*ND.DU2m/ND.TU2s/1000,'color',col,'LineWidth',1);
%
subplot(2,2,4);
for ii = 1 : num_MC
    plot(tspan*ND.TU2s/60,xhist(:,4,ii)*180/pi,'color',grey);
end
plot(tspan*ND.TU2s/60,x_nominal(:,4)*180/pi,'color',col,'LineWidth',1);

% Plot optimal control
figure(4); hold on; grid on;
xlabel('Time (min)','FontSize',15);
ylabel('Angle of Attack (deg)','FontSize',15);
set(gca,'FontSize',15);
stairs(tspan*ND.TU2s/60,u_nominal*180/pi,'color',col,'LineWidth',1);


% Compute errors between MC and nominal
error          = zeros(M, 4, num_MC);
final_error    = zeros(num_MC, 4);

for ii = 1 : num_MC
    
    error(:, :, ii)       = x_nominal - xhist(:, :, ii);
    final_error(ii, :)    = error(end, :, ii);
    
end

% Re-dimensionalize units
final_error(:, 1) = final_error(:, 1) * ND.DU2m / 1000;      % km
final_error(:, 2) = final_error(:, 2) * 180 / pi;            % deg
final_error(:, 3) = final_error(:, 3) * ND.DU2m / ND.TU2s;   % m/s
final_error(:, 4) = final_error(:, 4) * 180 / pi;            % deg

% Compute covariances and means
covariance_final   = cov(final_error);
mean_final         = mean(final_error)';
conf               = 0.9973;
Npl                = 100;
muf                = mean_final;
Sigmaf             = covariance_final;

% Ellipse semi major axes
length    = sqrt(-2 * log(1 - conf)) * sqrt(max(eig(covariance_final)));

% Plot styles
bl_scatter_style = {'o', 'MarkerEdgeColor', 0.4 * ones(3, 1)};
fdoc_scatter_style = {'k+'};
bl_plot_style = {'-.', 'Color', 0.25 * ones(3, 1)};
fdoc_plot_style = {'color',col};

% Plot final dispersions for baseline and FDOC
figure(5); hold on; grid on;
xlabel('Final altitude error (km)','FontSize',15);
ylabel('Final velocity error (m/s)','FontSize',15);
scatter(final_error(:,1),final_error(:,3),fdoc_scatter_style{:});
mean_alt_vel = [muf(1); muf(3)];
cov11 = Sigmaf(1,1);
cov13 = Sigmaf(1,3);
cov31 = Sigmaf(3,1);
cov33 = Sigmaf(3,3);
cov_alt_vel = [cov11 cov13; cov31 cov33];
conf_ellipse(mean_alt_vel,cov_alt_vel,Npl,conf,1,fdoc_plot_style{:});
set(gca,'FontSize',15);

figure(6); hold on; grid on;
xlabel('Final altitude error (km)','FontSize',15);
ylabel('Final longitude error (deg)','FontSize',15);
cov11 = Sigmaf(1,1);
cov12 = Sigmaf(1,2);
cov21 = Sigmaf(2,1);
cov22 = Sigmaf(2,2);
cov_alt_lon   = [cov11 cov12; cov21 cov22];
mean_alt_lon  = [muf(1); muf(2)];
scatter(final_error(:,1), final_error(:,2),fdoc_scatter_style{:});
conf_ellipse(mean_alt_lon,cov_alt_lon,Npl,conf,1,fdoc_plot_style{:});
set(gca,'FontSize',15);


end


