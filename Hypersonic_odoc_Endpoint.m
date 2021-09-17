function output = Hypersonic_odoc_Endpoint(input)

P     = input.auxdata.P;
Rm    = input.auxdata.Rm;
Qa    = input.auxdata.Qa;

rf  = input.phase.finalstate(1);

S1Hf   = input.phase.finalstate(5);
S2Hf   = input.phase.finalstate(6);
S3Hf   = input.phase.finalstate(7);
S4Hf   = input.phase.finalstate(8);

Sf = [S1Hf; S2Hf; S3Hf; S4Hf];

sensitivity_cost = 0.5*trace(Qa * Sf * P *Sf');

output.objective = -rf + sensitivity_cost;

end