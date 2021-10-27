function output = Hypersonic_Endpoint_v3_DOC(input)

alpha  = input.auxdata.alpha;
SigmaP = input.auxdata.SigmaP;
Qf     = input.auxdata.Qf;

Sfvec  = input.phase.finalstate(5 : 8); 
Sf     = reshape(Sfvec, 4, 1);

phif   = input.phase.finalstate(2);
        
sensitivity_cost = trace(Qf * Sf * SigmaP * Sf');

output.objective = -phif + alpha * sensitivity_cost;

end