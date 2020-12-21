function [Pn,Pms,Pmf] = Be_production_at_xypoint(y,z,Tshd)
% Calcuates scaling factos after Stone, 2000
% y - latitude in degree, z - elevation in m, Tshd, topographic shielding
% Richard Ott, 2019

% Spallation-induced Production Rate of 10Be (at/g/yr) assuming 'Lm" scaling framework (Borchers et al., 2015)
Pn_SLHL = 4.00; 

% Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011)
Pms_SLHL = 0.012; 

% Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011)
Pmf_SLHL = 0.039; 

% Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)to account for diminished cosmic ray attenuation with decreased atmospheric pressue 
pres=1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3.*z)))); 

%% scaling Stone 2000
% Scaling Equation Constants by Latitude (Stone,2000)
Lat = [0,10,20,30,40,50,60,90];
a = [31.8518,34.3699,40.3153,42.0983,56.7733,69.0720,71.8733,71.8733];
b = [250.3193,258.4759,308.9894,512.6857,649.1343,832.4566,863.1927,863.1927];
c = [-0.083393,-0.089807,-0.106248,-0.120551,-0.160859,-0.199252,-0.207069,-0.207069];
d = [7.4260e-5,7.9457e-5,9.4508e-5, 1.1752e-4,1.5463e-4,1.9391e-4,2.0127e-4,2.0127e-4];
e = [-2.2397e-8,-2.3697e-8,-2.8234e-8,-3.8809e-8,-5.0330e-8,-6.3653e-8,-6.6043e-8,-6.6043e-8];
m = [0.587,0.600,0.678,0.833,0.933,1.000,1.000,1.000];
    
A = interp1(Lat,a,y);
B = interp1(Lat,b,y);
C = interp1(Lat,c,y);
D = interp1(Lat,d,y);
E = interp1(Lat,e,y);
M = interp1(Lat,m,y);

%% Calculation 10Be production from spallation and muons

% 10Be production from Spallation, assuming neutron attenuation length in air of 150 g/cm2
Pn= Tshd.*Pn_SLHL.*(A+B.*exp(-pres/150)+C.*pres+D.*pres.^2+E.*pres.^3); % Stone, 2000

% 10Be production from Slow Muons, assuming (1) sea level pressure of 1013.25 mbar and 
%(2) muon attentuation length in air of 260 g/cm2 (Braucher et al., 2011)
Pms=Tshd.*Pms_SLHL.*exp((1013.25-pres)/260); 

% 10Be production from Fast Muons, assuming (1) sea level pressure of 1013.25 mbar and 
%(2) muon attentuation length in air of 510 g/cm2 (Braucher et al., 2011)
Pmf=Tshd.*Pmf_SLHL.*exp((1013.25-pres)/510); 

end