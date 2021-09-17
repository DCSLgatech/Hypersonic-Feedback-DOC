function output = Hypersonic_baseline_Endpoint(input)

rbarf  = input.phase.finalstate(1);

output.objective = -rbarf;

end