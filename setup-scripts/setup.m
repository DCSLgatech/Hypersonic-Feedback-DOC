%-------------------------------------------------------------------%
%--- Set all the constants, limits, initial and final conditions ---%
%-------------------------- appropriately --------------------------%
%-------------------------------------------------------------------%

function [C, IC, FC, LB, UB] = setup()

% C  - Constants
% IC - Initial conditions
% FC - Final conditions
% LB - Lower bounds
% UB - Upper bounds

%% Important constants

C.Re   = 6371000;           % Equatorial Radius of Earth (m)
C.S    = 150;               % Vehicle Reference Area (m^2)
C.H    = 7254.24;           % Density Scale Height (m)
C.rho0 = 1.225;             % Sea Level Atmospheric Density (kg/m^3)
C.g0   = 9.81;              % Gravity at sea level (m/s^2) 
C.mass = 38000;             % Vehicle Mass (kg)
C.kQ   = 9.4369 * 1e-5;     % heating rate constant

% Cost weights
C.Qf   = eye(4); %blkdiag(0, 1, 1, 0, 0, 0); % terminal output penalty
C.R    = 100;

% Path constraints

C.Qdotmax = 400000;   % Maximum heating rate allowed
C.qbarmax = 14500;    % Maximum dynamic pressure allowed
C.nmax    = 5 * C.g0; % Maximum normal load allowed/gravity

% Obtain Gravitational constant (m^3 / s^2)
C.mu   = C.g0 * C.Re^2;

%% Aerodynamic parameters

C.a0 = -0.20704;
C.a1 = 0.029244;
C.b0 = 0.07854;
C.b1 = -0.61592E-02;
C.b2 = 0.621408E-03;

%% Initial conditions

IC.time    = 0;                    % time (s)
IC.alt     = 121900;               % altitude (m) 121900
IC.speed   = 7626;                 % speed (m/s)
IC.phi     = -25 * pi / 180;       % longitude (rad)
IC.fpa     = -1.2493 * pi / 180;   % flight path angle (rad)
IC.rad     = IC.alt + C.Re; 
%IC.dr      = 0;
IC.Svec    = zeros(1, 4);

%% Final conditions

FC.alt     = 30480;             % altitude (m)
FC.speed   = 908.15;            % speed (m/s)
FC.phiMin  = -pi;               % minimum longitude (rad)
FC.phiMax  = pi;                % maximum longitude (rad)
FC.fpaMin  = -6 * pi/180;       % minimum flight path angle (rad)
FC.fpaMax  = 0 * pi/180;        % maximum flight path angle (rad)
FC.rad     = FC.alt + C.Re; 
FC.Pvec    = reshape(C.Qf, 1, 16);

%% Limits (or bounds)

% % Final time (s)
LB.tf  = 0; 
UB.tf  = 10000;

% radius (m)
LB.rad = C.Re; 
UB.rad = IC.rad;

% Longitude (rad)
LB.phi = -pi;         
UB.phi = pi;

% % Latitude (rad)
% LB.lat = -70 * pi / 180;  
% UB.lat = -LB.lat;

% Speed (m/s)
LB.speed = 10;        
UB.speed = 45000;

% Flight path angle (rad)
LB.fpa = -80 * pi / 180;  
UB.fpa = 80 * pi / 180;

% % Azimuth
% LB.azi = -180 * pi / 180; 
% UB.azi = -LB.azi;

% % Lift coefficient
% LB.CL  = -0.15;       
% UB.CL  = 0.8; 

% Angle of attack (rad)
LB.alpha = -80 * pi / 180;
UB.alpha = 80 * pi / 180;

% % Downrange distance (m)
% LB.dr = 0;
% UB.dr = 1E10;
% 
% % Bank angle (rad)
% LB.bank  = -90 * pi/180;
% UB.bank  = 90 * pi/180;

% % Control derivatives
% LB.CLdot   = -0.05; 
% UB.CLdot   = -LB.CLdot;
% 
% LB.bankdot = -5 * (pi / 180); 
% UB.bankdot = -LB.bankdot;

% Sensitivity Matrix
LB.Svec  = -1E5 * ones(1, 4); 
UB.Svec  = 1E5 * ones(1, 4);

% Riccati Matrix
LB.Pvec  = -1E5 * ones(1, 16);
UB.Pvec  = 1E5 * ones(1, 16);

end