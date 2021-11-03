function solution = obtain_solution_v3_FDOC(alpha, C, IC, FC, LB, UB, ND)

%% Auxiliary data for GPOPS-II

auxdata.alpha    = alpha;
auxdata.mubar    = C.mubar;
auxdata.rho0bar  = C.rho0bar;
auxdata.Hbar     = C.Hbar;
auxdata.Sbar     = C.Sbar;

auxdata.Re       = C.Re;
auxdata.g0       = C.g0;
auxdata.SigmaP   = C.SigmaP;
auxdata.R        = C.R;
auxdata.Qf       = C.Qf;
auxdata.Q        = C.Q;

auxdata.mass     = C.mass;
auxdata.kQ       = C.kQ;

auxdata.a0       = C.a0;
auxdata.a1       = C.a1;
auxdata.b0       = C.b0;
auxdata.b1       = C.b1;
auxdata.b2       = C.b2;

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
% Time bounds
bounds.phase.initialtime.lower = IC.tau;
bounds.phase.initialtime.upper = IC.tau;
bounds.phase.finaltime.lower   = LB.tauf;
bounds.phase.finaltime.upper   = UB.tauf;

% Initial state bounds
bounds.phase.initialstate.lower = [IC.rbar, IC.phi, IC.vbar, IC.fpa, IC.Svec, LB.Pvec];
bounds.phase.initialstate.upper = [IC.rbar, IC.phi, IC.vbar, IC.fpa, IC.Svec, UB.Pvec];

% State bounds
bounds.phase.state.lower = [LB.rbar, LB.phi, LB.vbar, LB.fpa, LB.Svec, LB.Pvec];                        
bounds.phase.state.upper = [UB.rbar, UB.phi, UB.vbar, UB.fpa, UB.Svec, UB.Pvec];
                        
% Final state bounds
bounds.phase.finalstate.lower = [FC.rbar, FC.phiMin, FC.vbar, FC.fpaMin, LB.Svec, FC.Pvec];
bounds.phase.finalstate.upper = [FC.rbar, FC.phiMax, FC.vbar, FC.fpaMax, UB.Svec, FC.Pvec];
                             
% Control bounds
bounds.phase.control.lower = LB.alpha;
bounds.phase.control.upper = UB.alpha;

% Bounds for constraint functions
bounds.phase.path.lower = zeros(1, 3);
bounds.phase.path.upper = [C.Qdotmax C.qbarmax C.nmax];

% Integral bounds
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 1E6;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%

% Component-wise guesses

% % % LOAD IN NOMINAL SOLUTION FOR BETTER PERFORMANCE % % %

nominal_loader  = load('fdoc_solution_a_10');
nominal         = nominal_loader.solution.phase;
nominal_time    = nominal.time;
nominal_state   = nominal.state;
nominal_control = nominal.control;


% tGuess               = [0; 1.44552];
% radGuess             = [IC.rbar; FC.rbar];
% phiGuess             = [IC.phi; 30 * pi / 180];
% speedGuess           = [IC.vbar; FC.vbar];
% fpaGuess             = [IC.fpa; -0.105];
% SvecGuess            = [zeros(1, 4); zeros(1,4)];
% PvecGuess            = [zeros(1, 16); FC.Pvec];

% Guesses for state, control and time
%guess.phase.state    = [radGuess, phiGuess, speedGuess, fpaGuess, SvecGuess, PvecGuess];
%guess.phase.control  = [0; 0];
%guess.phase.time     = tGuess;

guess.phase.state    = nominal_state;
guess.phase.control  = nominal_control;
guess.phase.time     = nominal_time;

guess.phase.integral = 1;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%

% if alpha == 0
%     mesh.method       = 'hp-LiuRao-Legendre';
%     % mesh.method       = 'hp-LiuRao';
%     % mesh.method       = 'hp-DarbyRao';
%     % mesh.method       = 'hp-PattersonRao';
%     mesh.maxiterations = 30;
%     mesh.colpointsmin  = 3;
%     mesh.colpointsmax  = 20;
%     mesh.tolerance     = 1E-03;
% else
    mesh.method       = 'hp-LiuRao-Legendre';
%     mesh.method       = 'hp-LiuRao';
%     mesh.method       = 'hp-DarbyRao';
%     mesh.method       = 'hp-PattersonRao';
    mesh.maxiterations = 30;
    mesh.colpointsmin  = 3;
    mesh.colpointsmax  = 20;
    mesh.tolerance     = 1e-3;
% end

% mesh.method       = 'hp-LiuRao-Legendre';
% % mesh.method       = 'hp-LiuRao';
% % mesh.method       = 'hp-DarbyRao';
% % mesh.method       = 'hp-PattersonRao';
% mesh.maxiterations = 30;
% mesh.colpointsmin  = 3;
% mesh.colpointsmax  = 20;
% mesh.tolerance     = 5e-5;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Hypersonic-Problem';
setup.functions.continuous           = @Hypersonic_Continuous_v3_FDOC;
setup.functions.endpoint             = @Hypersonic_Endpoint_v3_FDOC;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.derivatives.supplier           = 'sparseFD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output   = gpops2(setup);
solution = output.result.solution;

end