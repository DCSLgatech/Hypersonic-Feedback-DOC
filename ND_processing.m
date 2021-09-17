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
IC.rbar   = IC.rad * ND.m2DU; 
IC.vbar   = IC.speed * ND.m2DU / ND.s2TU;
IC.dbar   = IC.dran * ND.m2DU;

%% Final conditions

FC.vbarMin   = FC.speedMin * ND.m2DU / ND.s2TU;
FC.vbarMax   = FC.speedMax * ND.m2DU / ND.s2TU;
FC.dbarMin   = FC.dranMin * ND.m2DU; 
FC.dbarMax   = FC.dranMax * ND.m2DU; 

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
% Downrange (m)
LB.dbar = LB.dran * ND.m2DU; 
UB.dbar = UB.dran * ND.m2DU;

% Sensitivities
LB.sensibar = -1e8 * [1, 1, 1, 1]; 
UB.sensibar = -LB.sensibar;