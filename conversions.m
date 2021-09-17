% Useful conversions
% DU - Non-dimesnional distance unit
% TU - Non-dimesnional time unit
% m - meters
% s - seconds

ND.DU2m = C.Rm;
ND.m2DU = 1 / ND.DU2m;

ND.TU2s = sqrt(C.Rm^3 / C.mu);
ND.s2TU = 1 / ND.TU2s;