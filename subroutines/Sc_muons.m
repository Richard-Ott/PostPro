function [Scm] = Sc_muons(pressure)
% Returns the muonic sclaing factor need in many caluclations.
% muon attentuation length in air of 260 g/cm2 (Braucher et al., 2011)

pressure_SL = 1013.25;
lambda = 260;
Scm = (pressure_SL-pressure)./lambda; 
end

