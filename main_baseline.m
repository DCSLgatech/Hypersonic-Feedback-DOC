% ----------------------------------------------------------------------- %
% ------------- Hypersonic vehicle Trajectory Optimization -------------- %
% ------------------ Author: Joshua Pilipovsky  ------------------------- %
% ---------- (uncertain parameters: Scale height - H) ------------------- %
% ----------------------------------------------------------------------- %

close all;
clear all;
clc;

addpath('GPOPS-II-continuous');
addpath('GPOPS-II-endpoint');
addpath('GPOPS-II-obtain-solution');
addpath('linearizations');
addpath('post-processing');
addpath('setup-scripts');
addpath('solution-data-files');

% Choose alpha values for your desensitization 
% 0 implies no desensitization
alpha = 0;

% Load the set-up (constants, initial and final conditions, bounds)
[C, IC, FC, LB, UB] = setup();

% Obtain some important conversions
% DU - Non-dimensional distance unit
% TU - Non-dimesnional time unit
% m - meters
% s - seconds

ND.DU2m = C.Re;
ND.m2DU = 1 / ND.DU2m;

ND.TU2s = sqrt(C.Re / C.g0);
ND.s2TU = 1 / ND.TU2s;

% Nondimensionalization
% Non-dimensionalization
% bar represents non-dimendional quantities
% tau is the non-dimensional time

%% Important constants

C.Sbar    = C.S * ND.m2DU^2;         
C.Hbar    = C.H * ND.m2DU;           
C.rho0bar = C.rho0 / ND.m2DU^3;
C.mubar   = C.mu * ND.m2DU^3 / ND.s2TU^2;

%% Initial conditions

IC.tau    = IC.time * ND.s2TU;
%FC.tauf   = FC.time * ND.s2TU;

IC.rbar   = IC.rad * ND.m2DU; 
IC.vbar   = IC.speed * ND.m2DU / ND.s2TU;
% IC.dbar   = IC.dr * ND.m2DU;


%% Final conditions

if isfield(FC,'rad')
    FC.rbar  = FC.rad * ND.m2DU;
end
if isfield(FC,'speed')
    FC.vbar  = FC.speed * ND.m2DU / ND.s2TU;
end

%% Limits (or bounds)

% Final time (s)
LB.tauf  = LB.tf * ND.s2TU; 
UB.tauf  = UB.tf * ND.s2TU;

% radius (m)
LB.rbar = LB.rad * ND.m2DU; 
UB.rbar = UB.rad * ND.m2DU;

% Speed (m/s)
LB.vbar = LB.speed * ND.m2DU / ND.s2TU;        
UB.vbar = UB.speed * ND.m2DU / ND.s2TU;

% % Downrange (m)
% LB.dbar = LB.dr * ND.m2DU;
% UB.dbar = UB.dr * ND.m2DU;

% % Control derivatives
% LB.CLdotbar   = LB.CLdot / ND.s2TU; 
% UB.CLdotbar   = -LB.CLdotbar;
% 
% LB.bankdotbar = LB.bankdot / ND.s2TU; 
% UB.bankdotbar = -LB.bankdotbar;


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
tstart   = tic;
solution = obtain_solution_v3(alpha, C, IC, FC, LB, UB, ND);
simTime  = toc(tstart);

% Run MC analysis
[muf, Sigmaf] = MonteCarlo(solution,C,ND,IC,'k');

% Sensitivity analysis
loader = load('solution-data-files/open_loop_sensi');
sensi  = loader.S;
plotSensitivities(ND, solution.phase.time, sensi, 7, 'k');

