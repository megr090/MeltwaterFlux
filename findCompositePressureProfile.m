function [effpressure_composite] = findCompositePressureProfile(effpressure_outer,effpressure_inner)

effpressure_composite = effpressure_outer + effpressure_inner - effpressure_outer(1);

end

