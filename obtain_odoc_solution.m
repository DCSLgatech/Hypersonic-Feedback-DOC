function sol = obtain_odoc_solution(alpha, C, IC, FC, LB, UB)

%% Auxiliary data for GPOPS-II

auxdata.alpha    = alpha;
auxdata.mu       = C.mu;
auxdata.rho0     = C.rho0;
auxdata.H	     = C.H;
auxdata.S	     = C.S;

auxdata.CD       = C.CD;
auxdata.CL       = C.CL;
auxdata.Rm       = C.Rm;
auxdata.P        = C.P;
auxdata.Q        = C.Q;
auxdata.Qa       = C.Q * auxdata.alpha;

auxdata.mass     = C.mass;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
% Time bounds
bounds.phase.initialtime.lower = IC.time;
bounds.phase.initialtime.upper = IC.time;
bounds.phase.finaltime.lower   = LB.tf;
bounds.phase.finaltime.upper   = UB.tf;

% Initial state bounds
bounds.phase.initialstate.lower = [IC.rad,IC.speed,IC.fpa,IC.dran,zeros(1,4)];
bounds.phase.initialstate.upper = [IC.rad,IC.speed,IC.fpa,IC.dran,zeros(1,4)];

% State bounds
bounds.phase.state.lower = [LB.rad,LB.speed,LB.fpa,LB.dran,LB.sensi*ones(1,4)];                        
bounds.phase.state.upper = [UB.rad,UB.speed,UB.fpa,UB.dran,UB.sensi*ones(1,4)];
                        
% Final state bounds
bounds.phase.finalstate.lower = [LB.rad,FC.speedMin,FC.fpaMin,FC.dranMin,...
                                 LB.sensi*ones(1,4)];
bounds.phase.finalstate.upper = [UB.rad,FC.speedMax,FC.fpaMax,FC.dranMax,...
                                 UB.sensi*ones(1,4)];
                             
% Control bounds
bounds.phase.control.lower = [LB.cbank];
bounds.phase.control.upper = [UB.cbank];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%

% Component-wise guesses
tGuess              = [0; 216.5];
radGuess            = [IC.rad; 10530];
speedGuess          = [IC.speed; 348.8];
fpaGuess            = [IC.fpa; -23.5*pi/180];
dranGuess           = [IC.dran; 566000];
sensiGuess          = [zeros(1,4); -1.67, -0.05, -9.32*1e-5, -25.83];
cbankGuess          = [LB.cbank; UB.cbank];

% Guesses for state, control and time
guess.phase.state    = [radGuess, speedGuess,fpaGuess, dranGuess, sensiGuess];
guess.phase.control  = cbankGuess;
guess.phase.time     = tGuess;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%


mesh.method        = 'hp-LiuRao-Legendre';
%     mesh.method       = 'hp-LiuRao';
%     mesh.method       = 'hp-DarbyRao';
%     mesh.method       = 'hp-PattersonRao';
mesh.maxiterations = 30;
mesh.colpointsmin  = 3;
mesh.colpointsmax  = 40;
mesh.tolerance     = 1e-4;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Hypersonic-odoc-Problem';
setup.functions.continuous           = @Hypersonic_odoc_Continuous;
setup.functions.endpoint             = @Hypersonic_odoc_Endpoint;
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

%% Reconstruct solution using ode45
sol = reconstruct_solution(solution, C, ND, IC);

end