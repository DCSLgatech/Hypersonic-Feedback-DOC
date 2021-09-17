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

C.Rm   = 3397000;           % Equatorial Radius of Mars (m)
C.S    = 4.6^2 * pi/4;      % Vehicle Reference Area (m^2)
C.H    = 11000;             % Density Scale Height (m)
C.rho0 = 0.02;              % Sea Level Atmospheric Density (kg/m^3)
C.mass = 2920;              % Vehicle Mass (kg)

% Aerodynamic coefficients
C.CD   = C.mass/(115 * C.S); 
C.CLCD = 0.24;
C.CL   = C.CD * C.CLCD;

% Obtain Gravitational constant (m^3 / s^2)
C.mu   = 42828.29 * 1e9;

%% Initial conditions

IC.time  = 0;              % time (s)
IC.alt   = 125000;         % altitude (m) 121900
IC.speed = 5900;           % speed (m/s)
IC.fpa   = -15*pi/180; % flight path angle (rad)
IC.dran  = 0;              % Downrange (m)

% Converting altitude to radius
IC.rad   = IC.alt + C.Rm; 

%% Final conditions

FC.speedMin  = 0;               % minimum speed (m/s)
FC.speedMax  = 350;             % maximum speed (m/s)
FC.fpaMin    = -90*pi/180;       % minimum flight path angle (rad)
FC.fpaMax    = 90*pi/180;        % maximum flight path angle (rad)
FC.dranMin   = 100000;          % minimum downrange (m)
FC.dranMax   = 900000;          % maximum downrange (m)

%% Limits (or bounds)

% Final time (s)
LB.tf  = 0; 
UB.tf  = 3000;
% radius (m)
LB.rad = C.Rm; 
UB.rad = IC.rad+1e5;
% Speed (m/s)
LB.speed = 10;        
UB.speed = 45000;
% Flight path angle (rad)
LB.fpa = -80*pi/180;  
UB.fpa = -LB.fpa;
% Downrange (m)
LB.dran = 0;  
UB.dran = 1e6;

% Cosine bank angle
LB.cbank  = -0.8;
UB.cbank  = 0.8;

% Control derivatives
% LB.bankdot = -5*(pi/180); 
% UB.bankdot = -LB.bankdot;

% Sensitivities
LB.sensi = -1e8; 
UB.sensi = -LB.sensi;

LB.P  = -1E8 * ones(1, 10);
UB.P  = 1E8 * ones(1, 10);


% LB.P = [0,0,-30010.3119041599,0,0,300.733593946325,-4028454.95249168,0,...
%     -40910.3119041599,-4028454.95249168,0.72957795113294,0,0,0,0,...
%     0.69656219214308]; 
% UB.P = 1e12 * [486.900549428747,4222.24555549072,179770825.242065,132.429267483968,...
%     4222.24555549072,200000,131599846.748753,1793.88754366024,179770825.242065,...
%     131599846.748753,77067896855415.6,23004983.5668149,132.429267483968,...
%     1793.88754366024,23004983.5668149,50];

% h0 = 125 km
% v0 = 5800 m/s
% fpa = 15.5 deg
% 
% 
% dran = 600+ km

end