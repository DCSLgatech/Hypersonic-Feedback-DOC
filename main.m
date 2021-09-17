% ------------- Hypersonic vehicle Trajectory Optimization -------------- %
% ------------------ Author: Venkata Ramana Makkapati ------------------- %
% -- Run final time energy desensitization using Feedback DOC approach -- % 
% -------------- Considering longitudinal dynamics alone ---------------- %
% -------------- (uncertain parameters: Scale height - H) --------------- %

close all;
clear all;
clc;

% % add adigator to your path
% addpath(genpath('../adigator'))

% Load the set-up (constants, initial and final conditions, bounds)
[C, IC, FC, LB, UB] = setup();

% Obtain some important conversions
run conversions

% Nondimensionalization
run ND_processing

% Choose alpha values for your desensitization 
% Open loop DOC
odoc_alpha_values = [1e-3, 1e-1];  % Dont choose zero
% LQR DOC 
lqrdoc_alpha_values = [0.1];  % Dont choose zero

% Sample 100 parameters values within +-2% of nominal one for 
% Monte Carlo simulations
p_nom   = [C.Hbar];
n_trials = 100;

% Normal distribution
C.SigmaP  = (0.05 / 3 * C.Hbar).^2;
% C.Qf = 1e-12 * diag([1e2 * ND.DU2m^2, 2e8 * ND.DU2m^2 / ND.TU2s^2, 100 * 180 / pi, 5e4 * ND.DU2m^2]);
% C.R  = 1e2;
C.Qf = 1 * diag([10 10 10 10]);
C.R  = 100;

p_range = mvnrnd(p_nom, C.SigmaP, n_trials)';

% figure(6)
% hold on
% box on
% scatter(p_range * ND.DU2m/1000, zeros(1, n_trials), '+')
% title('Parameter Samples')
% xlabel('Scale height (km)')

%% Obtain baseline solution
baseline = obtain_baseline_solution(C, ND, IC, FC, LB, UB);

% plotting
figure(1)
subplot(2,2,1)
hold on
box on
plot(baseline.tplot, baseline.xplot(:,1), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('Altitude (km)', 'Interpreter', 'LaTeX')

subplot(2,2,2)
hold on
box on
plot(baseline.tplot, baseline.xplot(:,2), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('Speed (km/s)', 'Interpreter', 'LaTeX')

subplot(2,2,3)
hold on
box on
plot(baseline.tplot, baseline.xplot(:,3), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('Flight path angle (deg)', 'Interpreter', 'LaTeX')

subplot(2,2,4)
hold on
box on
plot(baseline.tplot, baseline.xplot(:,4), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('Downrange (km)', 'Interpreter', 'LaTeX')

figure(2)
hold on
box on
plot(baseline.tplot, baseline.cbank, 'k', 'Linewidth', 2)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$\cos \sigma$', 'Interpreter', 'LaTeX')

figure(3)
subplot(2,2,1)
hold on
box on
plot(baseline.tplot, abs(baseline.SHplot(:,1)), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$S_1$', 'Interpreter', 'LaTeX')

subplot(2,2,2)
hold on
box on
plot(baseline.tplot, abs(baseline.SHplot(:,2)), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$S_2$', 'Interpreter', 'LaTeX')

subplot(2,2,3)
hold on
box on
plot(baseline.tplot, abs(baseline.SHplot(:,3)), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$S_3$', 'Interpreter', 'LaTeX')

subplot(2,2,4)
hold on
box on
plot(baseline.tplot, abs(baseline.SHplot(:,4)), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$S_4$', 'Interpreter', 'LaTeX')

% Obtain initial guess for LQR DOC
baseline_lqr  = obtain_lqrgains(baseline, C);

baselinelqr_CLS = obtain_closedloopS(C, ND, baseline, baseline_lqr);

figure(3)
subplot(2,2,1)
plot(baselinelqr_CLS.t, abs(baselinelqr_CLS.SH(:,1)), 'k-.', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$S_1$', 'Interpreter', 'LaTeX')
legend('Open-loop', 'LQR-DOC Decoupled')

subplot(2,2,2)
plot(baselinelqr_CLS.t, abs(baselinelqr_CLS.SH(:,2)), 'k-.', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$S_2$', 'Interpreter', 'LaTeX')

subplot(2,2,3)
plot(baselinelqr_CLS.t, abs(baselinelqr_CLS.SH(:,3)), 'k-.', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$S_3$', 'Interpreter', 'LaTeX')

subplot(2,2,4)
plot(baselinelqr_CLS.t, abs(baselinelqr_CLS.SH(:,4)), 'k-.', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$S_4$', 'Interpreter', 'LaTeX')
sgtitle('Closed-loop sensitivities')

% figure(5)
% plot(baselinelqr_CLS.t, baselinelqr_CLS.CC, 'k', 'Linewidth', 1)
% xlabel('Time (s)', 'Interpreter', 'LaTeX')
% ylabel('Control covariance', 'Interpreter', 'LaTeX')


% % State and control covariance
% tcl = baselinelqr_CLS.t;
% PPcl = nan(4, 4, length(tcl));
% PPu = nan(length(tcl), 1);
% for k = 1:size(PPcl, 3)
%    Sk = baselinelqr_CLS.SH(k, :)';
%    PPcl(:, :, k) = Sk * C.SigmaP * Sk';
%    Kk = interp1(baseline_lqr.tau * ND.TU2s, baseline_lqr.K, tcl(k));
%    PPu(k) = Kk * PPcl(:, :, k) * Kk';
% end
% % Get nominal control on this time grid for use later
% urefcl = interp1(baseline.tplot, baseline.cbank, tcl);
% 
% % and open loop
% PPol = nan(4, 4, size(baseline.SHplot, 1));
% for k = 1:size(PPol, 3)
%    Sk = baseline.SHplot(k, :)';
%    PPol(:, :, k) = Sk * C.SigmaP * Sk';
% end
% tol = baseline.tplot;
% 
% % Plot state covariance
% scales = [1e-3, 1e-3, 180 / pi, 1e-3, 1];
% labels = {'Altitude error (km)', 'Velocity error (km/s)', ...
%    'FPA error (deg)', 'Range error (km)', 'Cosine bank'};
% 
% figure
% for i = 1:4
%    subplot(2, 2, i)
%    hold on
%    % Open loop
%    erri = scales(i) * 3 * sqrt(squeeze(PPol(i, i, :)));
%    plol = plot(tol, +erri, 'k');
%    plot(tol, -erri, 'k')
% %     axis manual
%    % Closed loop
%    erri = scales(i) * 3 * sqrt(squeeze(PPcl(i, i, :)));
%    plcl = plot(tcl, +erri, 'k-.');
%    plot(tcl, -erri, 'k-.')
%    if i == 1
%        legend([plol, plcl], 'Open-loop', 'Closed-loop')
%    end
% end
% for i = 1:4
%    subplot(2, 2, i)
%    grid on
%    box on
%    xlabel('Time (s)')
%    ylabel(labels{i})
% end
% 
% figure
% hold on
% err = scales(5) * 3 * sqrt(PPu);
% plot(tcl, urefcl, 'k')
% plot(tcl, urefcl + err, 'k--')
% plot(tcl, urefcl - err, 'k--')
% ylim([-1.5, 1.5])


%% LQR DOC
% Sensitivity penalty for different alpha values
lqrdoc_sensi_cost = zeros(numel(lqrdoc_alpha_values),1);

% Run desensitization for different values of alpha
for i = 1:numel(lqrdoc_alpha_values)
    
    alpha = lqrdoc_alpha_values(i);
    
    % Obtain solution from GPOPS-II
    lqrdoc_solution(i) = obtain_lqrdoc_solution(alpha, C, ND, IC, FC, LB, UB);
end 


%% Open-loop DOC
% Sensitivity penalty for different alpha values
odoc_sensi_cost = zeros(numel(odoc_alpha_values),1);

colormap = 'brgmcy';

% Run desensitization for different values of alpha
for i = 1:numel(odoc_alpha_values)
    
    alpha = odoc_alpha_values(i);
    
    % Obtain solution from GPOPS-II
    odoc_solution(i) = obtain_odoc_solution(alpha, C, IC, FC, LB, UB);
    
    % plotting
    figure(1)
    subplot(2,2,1)
    plot(odoc_solution(i).t, odoc_solution(i).xplot(:,1), '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('Altitude (km)', 'Interpreter', 'LaTeX')
    
    subplot(2,2,2)
    plot(odoc_solution(i).t, odoc_solution(i).xplot(:,2), '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('Speed (km/s)', 'Interpreter', 'LaTeX')
    
    subplot(2,2,3)
    plot(odoc_solution(i).t, odoc_solution(i).xplot(:,3), '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('Flight path angle (deg)', 'Interpreter', 'LaTeX')
    
    subplot(2,2,4)
    plot(odoc_solution(i).t, odoc_solution(i).xplot(:,4), '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('Downrange (km)', 'Interpreter', 'LaTeX')
    
    figure(2)
    plot(odoc_solution(i).t, odoc_solution(i).cbank, '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$\cos \sigma$', 'Interpreter', 'LaTeX')
    
    figure(3)
    subplot(2,2,1)
    plot(odoc_solution(i).t, odoc_solution(i).SH(:,1), '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$S_1$', 'Interpreter', 'LaTeX')
    
    subplot(2,2,2)
    plot(odoc_solution(i).t, odoc_solution(i).SH(:,2), '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$S_2$', 'Interpreter', 'LaTeX')
    
    subplot(2,2,3)
    plot(odoc_solution(i).t, odoc_solution(i).SH(:,3), '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$S_3$', 'Interpreter', 'LaTeX')
    
    subplot(2,2,4)
    plot(odoc_solution(i).t, odoc_solution(i).SH(:,4), '--', 'color', colormap(i), 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$S_4$', 'Interpreter', 'LaTeX')
end