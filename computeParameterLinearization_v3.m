function D = computeParameterLinearization_v3(x, u, p, UQ)

% State, control, and parameters
vbar  = x(3);
Hbar  = p(1);

% Useful quantities
Dbar   = UQ(1);
Lbar   = UQ(2);
mass   = UQ(3);
hbar   = UQ(9);

% Compute analytical Jacobian
df1dH   = 0;
df2dH   = 0;
df3dH   = -Dbar * hbar / (mass * Hbar ^ 2);
df4dH   = Lbar * hbar / (mass * vbar * Hbar ^ 2);
dfdH    = [df1dH; df2dH; df3dH; df4dH];

D       = dfdH;
     
end