% ----------------------------------------------------------------------- %
% ------------- Hypersonic vehicle Trajectory Optimization -------------- %
% ------------------ Author: Joshua Pilipovsky  ------------------------- %
% ---------- (uncertain parameters: Scale height - H) ------------------- %
% ----------------------------------------------------------------------- %

close all;
clear all;
clc;

% Choose alpha values for your desensitization 
% 0 implies no desensitization
alpha = 0;

% Load the set-up (constants, initial and final conditions, bounds)
[C, IC, FC, LB, UB] = setup();

% Obtain some important conversions
run conversions

% Nondimensionalization
run ND_processing

% Sample 100 parameters values within +-2% of nominal one for 
% Monte Carlo simulations
%C.SigmaP = (0.02 / 3*[C.Hbar, 0;0, C.CD0]).^2;
C.SigmaP = (0.05 / 3)^2 * C.Hbar^2;
%p_nom    = [C.Hbar; C.CD0];
p_nom    = C.Hbar;
n_trials = 1000;

% Normal distribution
p_range = mvnrnd(p_nom, C.SigmaP, n_trials)';
    
% Obtain solution from GPOPS-II
solution = obtain_solution_v3(alpha, C, IC, FC, LB, UB, ND);

% Run MC analysis
[muf, Sigmaf] = MonteCarlo(solution,C,ND,IC,'k');

% Sensitivity analysis
loader = load('open_loop_sensi');
sensi  = loader.S;
plotSensitivities(solution.phase.time, sensi, 7, 'k');

