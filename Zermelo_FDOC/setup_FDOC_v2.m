%-------------------------------------------------------------------%
%--- Set all the constants, limits, initial and final conditions ---%
%-------------------------- appropriately --------------------------%
%-------------------------------------------------------------------%

function [C, IC, FC, LB, UB] = setup_FDOC_v2()

% C  - Constants
% IC - Initial conditions
% FC - Final conditions
% LB - Lower bounds
% UB - Upper bounds

%% Important constants

C.p1      = 10;
C.p2      = -0.5;
C.p       = [C.p1; C.p2];
eta1      = 0.1;
eta2      = 0.05;
sigma1    = eta1 / 3 * C.p1;
sigma2    = eta2 / 3 * C.p2;
C.SigmaP  = blkdiag(sigma1^2, sigma2^2);
C.Qf      = blkdiag(10, 100);
C.Q       = eye(2);
C.R       = 20;
C.alpha   = 50;

%% Initial conditions

IC.time  = 0;              % time
IC.x1    = 0;              % downstream 
IC.x2    = 0;              % upstream 
IC.Svec  = zeros(1,4);     % sensitivity matrix   (vectorized)

%% Final conditions

FC.x2    = 0;
FC.Pvec  = reshape(C.Qf,1,4);

%% Limits (or bounds)

% Final time
LB.tf     = 1; 
UB.tf     = 1;

% Downstream
LB.x1     = 0;
UB.x1     = 10;

% Upstream
LB.x2     = 0;
UB.x2     = 5;

% Steering angle
LB.u      = -pi/2;
UB.u      = pi/2;

% Sensitivity Matrix
LB.Svec  = -1E5 * ones(1,4); 
UB.Svec  = 1E5 * ones(1,4);

% Riccati Matrix
LB.Pvec  = -1E5 * ones(1,4);
UB.Pvec  = 1E5 * ones(1,4);

end