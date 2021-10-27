% Useful conversions
% DU - Non-dimesnional distance unit
% TU - Non-dimesnional time unit
% m - meters
% s - seconds

ND.DU2m = C.Re;
ND.m2DU = 1 / ND.DU2m;

ND.TU2s = sqrt(C.Re / C.g0);
ND.s2TU = 1 / ND.TU2s;