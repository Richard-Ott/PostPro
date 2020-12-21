function Leff = neutron_att_length(pressure,Z,lat,lon)
% This script calculates the neutron attenuation length according to CRONUS
% 2.0, Marrero (2016) and Sato (2008) with a dependency on time and 
% air-pressure.
% Magnetic time evolution parameters for neutron attenuation length 
% calculation are from CRONUS 2.0 (Marrero, 2016).
% Input:        -
% Richard Ott, 2020

addpath '.\subroutines\CRONUS_subroutines'

%% NEUTRON ATTENUATION LENGTH ---------------------------------------------
load pmag_consts.mat                  % get magentic field constants from CRONUS
% prepare sample file for rigiditycutoff calculation
sample.lat = lat;
sample.long = lon;
sample.scaling = 'sa'; 
sample.elevation = Z;
sample.pressure  = pressure;

% rigidity cutoff
tdsf = get_tdsf(sample,pmag_consts);
Rc = mean(tdsf.Rc_Sa);

Leff = rawattenuationlength(pressure,Rc);
end

