function B = computeControlLinearization_v3(x, u, p, UQ)

% State, control, and parameters
vbar  = x(3);
alpha = u;

% Useful quantities
mass   = UQ(3);
qbar   = UQ(5);
a1     = UQ(6);
b1     = UQ(7);
b2     = UQ(8);
Sbar   = UQ(10);

% Compute gradients of drag and lift
dD_dalpha = (180 / pi) * qbar * Sbar * (b1 + 2 * b2 * 180 / pi * alpha);
dL_dalpha = (180 / pi) * qbar * Sbar * a1;

% Compute analytical Jacobian
df1du   = 0;
df2du   = 0;
df3du   = -(1 / mass) * dD_dalpha;
df4du   = (1 / (mass * vbar)) * dL_dalpha;
dfdu    = [df1du; df2du; df3du; df4du];

B       = dfdu;
     
end