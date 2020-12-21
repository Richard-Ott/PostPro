function out = DepthProd_Heisinger_Balco2017(z,lat,lon,alt,nuc,att_method)
%%%%%%%%%%%%%%%%
% Function to calculate the 10Be production rate at depth
% Alt. and lat. scaling according to Lal1991/Stone2000
% Uses:	 - an exponential profile for neutrons and 
%		 - a fast and stopping muons production according to Heisinger's formulation. This part is directly adapted from Balco's and CRONUS P_mu_total.m file
%			 includes updated parameter values from Balco et al. 2017 for "Model 1A" with alpha = 1
%%%%%%%%%%%%%%%% 23.06.2017


%INPUT -----------------------------------------------------------------------------------------------------
%	z  						%postive vector of depth at which the production rate should be determined (g/cm2)
%	lat						%latitude of the site for Lal1991/Stone2000 scaling (Â°)
%   alt						%altitude of the site for Lal1991/Stone2000 scaling (m. asl)
%	nuc						%nuclide: 10, 26 or 14 
%   att_method              'manual' for manually adding the attenuation
%                           length or 'calc' for CRONUS based calculation
%                           of effective attenuation length (Marrero, 2016)
%OUTPUT ----------------------------------------------------------------------------------------------------
%	out	matrix of production rate for each depth in input z and each pathway (at/g/yr)
%Author: M. Lupker, migrated to Matlab by Richard Ott, 2020

%CONSTANTS -------------------------------------------------------------------------------------------------
if nuc  == 10
	Pn_SLHL = 4.01;			%neutron spallation SLHL production rate (10Be = 4.01 for Lal1991/Stone2000 scaling - according to Phillips et al., 2016) (at/g/yr)
	k_neg = 0.0002458;		%summary probability of negative muon capture (atoms/muon) (Heisinger et al., 2002 formulation: kneg = fc x fd x f*  with f* = 0.00191 according to Balco, 2017)
	sigma190 = 53.2e-30;	%190 Gev cross-section for fast muon production (cm2) (Blaco, 2007 = 53.2 Âµb / 1 Âµb = 1e-30cm2)
% 	att_neut = 160; 		%neutron attenuation length (g/cm2)
	Natoms = 2.006e22;		%atom number density of target atom (atoms/G)
	aalpha = 1;				%exponent for the energy dependence of the cross-section
	
elseif nuc  == 26
	Pn_SLHL = 27.93;		%neutron spallation SLHL production rate (26Al = 27.93 for Lal1991/Stone2000 scaling - according to Borchers et al., 2016) (at/g/yr)
	k_neg = 0.0002582;		%summary probability of negative muon capture (atoms/muon) (Heisinger et al., 2002 formulation: kneg = fc x fd x f*  with f* = 0.00133 according to Balco, 2017)
	sigma190 = 739e-30;		%190 Gev cross-section for fast muon production (cm2) (Blaco, 2017 = 739 Âµb / 1 Âµb = 1e-30cm2)
% 	att_neut = 160; 		%neutron attenuation length (g/cm2)
	Natoms = 1.003e22;		%atom number density of target atom (atoms/G)
	aalpha = 1;				%exponent for the energy dependence of the cross-section
elseif nuc  == 14
	Pn_SLHL = 12.24;		%neutron spallation SLHL production rate (26Al = 27.93 for Lal1991/Stone2000 scaling - according to Borchers et al., 2016) (at/g/yr)
	k_neg = 0.0149282;		%summary probability of negative muon capture (atoms/muon) (Heisinger et al., 2002 formulation: kneg = fc x fd x f*  with f* = 0.116 according to Lupker et al., 2015 & Balco, 2017)
	sigma190 = 450e-30;		%190 Gev cross-section for fast muon production (cm2) (Heisinger et al., 2002 = 450 Âµb / 1 Âµb = 1e-30cm2)
% 	att_neut = 160; 		%neutron attenuation length (g/cm2)
	Natoms = 2.006e22;		%atom number density of target atom (atoms/G)
	aalpha = 1;				%exponent for the energy dependence of the cross-section	
else
    disp("Nuclide not recognized")
end


if strcmpi(att_method,'manual')
    att_neut = input('Please enter your effective attenuation length for neutrons in g/cm� ');
elseif strcmpi(att_method,'calc')
    % Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)to account for diminished cosmic ray attenuation with decreased atmospheric pressue 
    pressure = 1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3.*alt)))); 
    att_neut = neutron_att_length(pressure,alt,lat,lon);
else
    error('read the instrcutions, man...')
end


%LAL1991/STONE2000 SCALING factors -------------------------------------------------------------------------

Latitude = [0,10,20,30,40,50,60];			%Scaling coefficients as a function of latitude (Stone et al, 2000)
A		 = [31.8518,34.3699,40.3153,42.0983,56.7733,69.0720,71.8733];
B		= [250.3193,258.4759,308.9894,512.6857,649.1343,832.4566,863.1927];
C		= [-0.083393,-0.089807,-0.106248,-0.120551,-0.160859,-0.199252,-0.207069];
D		= [7.4260e-5,7.9457e-5,9.4508e-5, 1.1752e-4,1.5463e-4,1.9391e-4,2.0127e-4];
E		= [-2.2397e-8,-2.3697e-8,-2.8234e-8,-3.8809e-8,-5.0330e-8,-6.3653e-8,-6.6043e-8];
M		= [0.587,0.600,0.678,0.833,0.933,1.000,1.000];

PRES=1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*alt)))); %Mean atmospheric pressure

a = interp1(Latitude,A,lat);  	%Interpolation of Stones empirical parameters according to input lat
b = interp1(Latitude,B,lat);
c = interp1(Latitude,C,lat);
d = interp1(Latitude,D,lat);
e = interp1(Latitude,E,lat);
m = interp1(Latitude,M,lat);

Fn_St = (a+b*exp(-PRES/150)+c*PRES+d*PRES^2+e*PRES^3);
Fm_St = m*exp((1013.25-PRES)/242);

%NEUTRON PRODUCTION RATE -----------------------------------------------------------------------------------
P_neut = Pn_SLHL * Fn_St * exp(-z./att_neut);

%MUON PRODUCTION -------------------------------------------------------------------------------------------
% Formulation after Heisinger's papers and directly adapted from Balco & CRONUS P_mu_total.m file

% Main function

% figure the atmospheric depth in g/cm2																						
h = 1013.25 * exp(((-0.03417) / 6.5e-3) * (log(288.15) - log(288.15 - (6.5e-3 * alt))));			%local atmospheric pressure according to global model
H = (1013.25 - h) * 1.019716;																		%local atmospheric depth

% vertical flux at SLHL
a = 258.5 * (100 ^ 2.66);
b = 75 * (100 ^ 1.66);
phi_vert_slhl = (a ./ ((z + 21000) .* (((z + 1000) .^ 1.66) + b))) .* exp(-5.5e-6 .* z);

% find the stopping rate of vertical muons at SLHL
R_vert_slhl = Rv0(z);

% find the stopping rate of vertical muons at site
R_vert_site = R_vert_slhl .* exp(H ./ LZ(z));

%flux of vertical muons at site
phi_vert_site = nan(1,length(z));

phi_vert_func = @(x) Rv0(x) .* exp(H ./ LZ(x));  
for i= 1:length(z)
	%integrate
	%ends at 200'001 g/cm2 to avoid being asked for an zero range integration
	%integration relative tolerance 1 part per 10^4
	tol = phi_vert_slhl(i) * 1e-4;
	phi_vert_site(i) = integral(phi_vert_func,z(i),(2e5 + 1),'AbsTol', tol);	
end

%invariant flux at 2e5 g/cm2 depth - constant of integration
phi_200k = (a / ((2e5 + 21000) * (((2e5 + 1000) ^ 1.66) + b))) * exp(-5.5e-6 * 2e5);
phi_vert_site = phi_vert_site + phi_200k;

%angular distribution exponent
nofz = 3.21 - 0.297 .* log((z + H) / 100 + 42) + 1.21e-5 .* (z + H);
%derivative
dndz = (-0.297 / 100) ./ ((z + H) ./ 100 + 42) + 1.21e-5;
phi_temp = phi_vert_site * 2 * pi ./ (nofz + 1);

%convert from muons/cm2/s to muons/cm2/yr
phi = phi_temp * 60 * 60 * 24 * 365;

%find the total stopping rate of muons at site
R_temp = (2 * pi ./ (nofz + 1)) .* R_vert_site - phi_vert_site .* (-2 * pi .* ((nofz + 1).^-2)) .* dndz;

%convert from muons/g/s to negative muons/g/yr
R = R_temp .* 0.44 .* 60 .* 60 .* 24 .* 365;

%Compute production rates
%Beta = 0.846 - 0.015 * log((z/100)+1) + 0.003139 * (log((z/100)+1)^2)
Beta = 1; 																				%beta = 1 for alpha = 1
Ebar = 7.6 + 321.7 .* (1 - exp(-8.059e-6 .* z)) + 50.7 .* (1-exp(-5.05e-7 .* z));

%internally defined constants

sigma0 = sigma190 / (190 ^ aalpha);

%fast muon production
P_fast = phi .* Beta .* (Ebar .^ aalpha) .* sigma0 .* Natoms;

%negative muon capture
P_neg = R .* k_neg;

%output the 10Be production of each pathways (at/g/yr)

out = [P_neut',P_neg',P_fast'];
end


%Auxiliary sub-function
function out = LZ(z)		%effective atmospheric attenuation length for muons of range Z (data table from Groom and others 2001)
		%units are range in g/cm2 (column 2) and momentum MeV/c (column 1)
		data =                  [4.704e1,8.516e-1;
    							5.616e1,1.542e0;
    							6.802e1,2.866e0;
    							8.509e1,5.698e0;
    							1.003e2,9.145e0;
    							1.527e2,2.676e1;
    							1.764e2,3.696e1;
	    						2.218e2,5.879e1;
   		 					    2.868e2,9.332e1;
    							3.917e2,1.524e2;
    							4.945e2,2.115e2;
    							8.995e2,4.418e2;
    							1.101e3,5.534e2;
    							1.502e3,7.712e2;
    							2.103e3,1.088e3;
    							3.104e3,1.599e3;
    							4.104e3,2.095e3;
    							8.105e3,3.998e3;
    							1.011e4,4.920e3;
    							1.411e4,6.724e3;
    							2.011e4,9.360e3;
    							3.011e4,1.362e4;
    							4.011e4,1.776e4;
    							8.011e4,3.343e4;
    							1.001e5,4.084e4;
    							1.401e5,5.495e4;
    							2.001e5,7.459e4;
    							3.001e5,1.040e5;
    							4.001e5,1.302e5;
    							8.001e5,2.129e5];
	
	%deal with zero depth situation
	too_low = find(z < 1);
	z(too_low) = ones(1,length(too_low));

	%obtain momenta with log-linear interpolation
	P_MeVc = exp(interp1(log(data(:,2)),log(data(:,1)),log(z)));

	out = (263 + 150 * P_MeVc / 1000);
end

%Auxiliary sub-function
function out = Rv0(z)		%stopping rate of vertically traveling muons as a function of depth z at SLHL
		a = exp(-5.5e-6 .* z);
		b = z + 21000;
		c = (z + 1000).^ 1.66 + 1.567e5;
		dadz = -5.5e-6 * exp(-5.5e-6 .* z);
		dbdz = 1;
		dcdz = 1.66 .* (z + 1000).^ 0.66;
		out = (-5.401e7 .* (b .* c .* dadz - a .* (c .* dbdz+b .* dcdz))./(b.^2 .* c.^2));
end