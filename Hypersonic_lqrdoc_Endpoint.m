function output = Hypersonic_lqrdoc_Endpoint(input)

alpha  = input.auxdata.alpha;
SigmaP = input.auxdata.SigmaP;
Qf     = input.auxdata.Qf;


rbarf  = input.phase.finalstate(1);

S1Hf   = input.phase.finalstate(5);
S2Hf   = input.phase.finalstate(6);
S3Hf   = input.phase.finalstate(7);
S4Hf   = input.phase.finalstate(8);

Sf = [S1Hf; S2Hf; S3Hf; S4Hf];

sensitivity_cost = 0.5*trace(Qf * Sf * SigmaP * Sf');

output.objective = -rbarf + alpha * (sensitivity_cost + input.phase.integral);

end