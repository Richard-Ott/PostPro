function Leff = neutron_att_length(lat,lon,elevation)
% This script calculates the neutron attenuation length according to CRONUS
% 2.0, Marrero (2016) and Sato (2008) with a dependency on time and 
% air-pressure.
% Magnetic time evolution parameters for neutron attenuation length 
% calculation are from CRONUS 2.0 (Marrero, 2016).
% Input:        - lat,lon, elevation
% Richard Ott, 2021

% addpath '.\CRONUS_subroutines'

%% NEUTRON ATTENUATION LENGTH ---------------------------------------------

% Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)to account for diminished cosmic ray attenuation with decreased atmospheric pressue 
pressure = 1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*elevation)))); 


Leff = nan(length(pressure),1);

load pmag_consts.mat                  % get magentic field constants from CRONUS
% prepare sample file for rigiditycutoff calculation
sample.lat = lat;
sample.long = lon; 
sample.scaling = 'sa'; 
sample.elevation = elevation;
sample.pressure = pressure;

% rigidity cutoff
tdsf = get_tdsf(sample,pmag_consts);
Rc = mean(tdsf.Rc_Sa);

Leff = rawattenuationlength(single(pressure),single(Rc)); 

end

