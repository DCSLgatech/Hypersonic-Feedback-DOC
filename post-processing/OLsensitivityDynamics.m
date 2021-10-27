function Sdot = OLsensitivityDynamics(S, A_hist, D_hist, t, ts)

A = interp1(ts, A_hist, t);
A = [A(:, :, 1); A(:, :, 2); A(:, :, 3); A(:, :, 4)];
D = interp1(ts, D_hist, t)';

Sdot = A * S + D;

end