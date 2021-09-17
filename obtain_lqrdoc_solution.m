function sol = obtain_lqrdoc_solution(alpha, C, ND, IC, FC, LB, UB)

%% Auxiliary data for GPOPS-II

auxdata.alpha    = alpha;
auxdata.mubar    = C.mubar;
auxdata.rho0bar  = C.rho0bar;
auxdata.Hbar     = C.Hbar;
auxdata.Sbar     = C.Sbar;

auxdata.CD       = C.CD;
auxdata.CL       = C.CL;
auxdata.Rm       = C.Rm;
auxdata.SigmaP   = C.SigmaP;
auxdata.R        = C.R;
auxdata.Qf       = C.Qf;

auxdata.mass     = C.mass;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
% Time bounds
bounds.phase.initialtime.lower = IC.tau;
bounds.phase.initialtime.upper = IC.tau;
bounds.phase.finaltime.lower   = LB.tauf;
bounds.phase.finaltime.upper   = UB.tauf;

% Initial state bounds
bounds.phase.initialstate.lower = [IC.rbar,IC.vbar,IC.fpa,IC.dbar,zeros(1,4),LB.P];
bounds.phase.initialstate.upper = [IC.rbar,IC.vbar,IC.fpa,IC.dbar,zeros(1,4),UB.P];

% State bounds
bounds.phase.state.lower = [LB.rbar,LB.vbar,LB.fpa,LB.dbar,LB.sensibar,LB.P];                        
bounds.phase.state.upper = [UB.rbar,UB.vbar,UB.fpa,UB.dbar,UB.sensibar,UB.P];
                        
% Final state bounds
bounds.phase.finalstate.lower = [LB.rbar,FC.vbarMin,FC.fpaMin,FC.dbarMin,...
                                 LB.sensibar, [0 0 0 0 1 0 0 0 0 0]];
bounds.phase.finalstate.upper = [UB.rbar,FC.vbarMax,FC.fpaMax,FC.dbarMax,...
                                 UB.sensibar, [0 0 0 0 1 0 0 0 0 0]];
                             
% Control bounds
bounds.phase.control.lower = [LB.cbank];
bounds.phase.control.upper = [UB.cbank];

bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 1e6;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%

% Component-wise guesses
tGuess              = [0; 216.5*ND.s2TU];
radGuess            = [IC.rbar; (C.Rm + 10530)*ND.m2DU];
speedGuess          = [IC.vbar; 348.8*ND.m2DU/ND.s2TU];
fpaGuess            = [IC.fpa; -23.5*pi/180];
dranGuess           = [IC.dran; 566000*ND.m2DU];
% sensiGuess          = [zeros(1,4); -1.67, -0.05/ND.s2TU, -9.32*1e-5/ND.m2DU, -25.83];
% PGuess              = [8.0655, -0.0591, 1.1315, 0, 0.0004, -0.0083, 0, 0.1587, 0, 0; 0 0 0 0 1 0 0 0 0 0];
sensiGuess          = [zeros(1,4);-1.0293  -38.0183 -258.0774  -23.8445];
PGuess              = [255.6465   -1.8049   35.7895   -0.0000    0.0127   -0.2527    0.0000    5.0104    0.0000   -0.0000; 0 0 0 0 1 0 0 0 0 0];
cbankGuess          = [-0.8; 0.8];

% Guesses for state, control and time
guess.phase.state    = [radGuess, speedGuess,fpaGuess, dranGuess, sensiGuess, PGuess];
guess.phase.control  = cbankGuess;
guess.phase.time     = tGuess;
guess.phase.integral = 1;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%


mesh.method       = 'hp-LiuRao-Legendre';
%     mesh.method       = 'hp-LiuRao';
%     mesh.method       = 'hp-DarbyRao';
%     mesh.method       = 'hp-PattersonRao';
mesh.maxiterations = 30;
mesh.colpointsmin  = 3;
mesh.colpointsmax  = 20;
mesh.tolerance     = 1e-3;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Hypersonic-lqrdoc-Problem';
setup.functions.continuous           = @Hypersonic_lqrdoc_Continuous;
setup.functions.endpoint             = @Hypersonic_lqrdoc_Endpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output   = gpops2(setup);
solution = output.result.solution;

plotLQRDOC(solution, ND, C);

%% Reconstruct solution using ode45
sol = reconstruct_solution(solution, C, ND, IC);

end