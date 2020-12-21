% This code calculates depth profiles of post burial production for 10Be
% and 36Cl (could easily expanded to work for 26Al and 14C). The
% calculations are based on geochronologic anchor ages. For the sites with
% geochronologic ages the site specific parameters (latitude, elevation,
% shielding, if applicable chemistry etc) need to be specified. The code
% then takes the one or more samples and runs Monte Carlo simulation of
% deposition, where deposition of a random amount of layers with certain
% thickness and time intervals follow power law distributions, as suggested
% by Gomez et al. 2002.
% Calculates the attenuation length for spallation following CRONUS and
% Sato (2008) with a dependency on time and air-pressure.
% This code still uses exponentials for muons production. I'm having
% problems understanding how to modify the Balco 2017 code to get
% attenuation lengths for fast and slow muons.
% The code is based on an R-script by M. Lupker.
%
% Richard Ott, 2020
clc
clear
close all

addpath '.\subroutines'
addpath '..\data'

% USER CHOICE ----------------------------------------------------------- %
nuclide = '10Be';       % Choose '10Be' or '36Cl'
export = 0;             % do you want to save figures and model.data
n = 5e3;                % number of runs

%% assign data and constants -------------------------------------------- %
p_thickness = 2.06; 	% parameter for power law distribution of sediment layer thickness (Gomez et al., 2002)
p_time = 1.4;		    % parameter for power law distribution of sediment layer waiting time (Gomez et al., 2002)
e = 0;                  % erosion rate (mm/yr)
att_neut_uncert = 0.05;	% relative uncertainty on the attenuation length for neutrons
att_sm_uncert = 0.5;	% relative uncertainty on the attenuation length for stopping muons
att_fm_uncert = 0.5;    % relative uncertainty on the attnunation length for fast muons
switch nuclide
    case '10Be'
        addpath '..\..\..\..\Crete\Cretan_fans\data'
        [data,~,raw] = xlsread('Samples_09-19.xlsx','10Be'); %load sample data
        ind = input('Enter the number of the sample you want to run (the row within the data file) ');
        depth = data(ind,11)*100;        % depth of sample (cm)
        depth_uncert = data(ind,12)*100; % uncertainty in overburden (cm)
        age = data(ind,13)*1e3;         % age of sample   (yrs)
        age_uncert = data(ind,14)*1e3;  % age uncertainty (yrs)
        lat = data(ind,3);              % latitude in degree, must be same for all depths
        lon = data(ind,4);              % longitude in degree (needed for cutoff rigidity)
        elev = data(ind,5);             % elevation in m
        Tshd = data(ind,15);            % topographic shielding factor
        name = raw(1+ind,1);

        % CONSTANTS ----------------------------------------------------- %
        rho = 2.2;              % density g/cm³
        lambda = log(2)/1.39e6;	% decay constant for 10Be
        air_pressure = 1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*elev)))); 
        
        att_neut = neutron_att_length(air_pressure,elev,lat,lon); % following CRONUS 2.0, Marrero 2016
        att_sm = 1500;
        att_fm = 4320;
        % UNCERTAINTIES ------------------------------------------------- %
        rho_uncert = 0.1;	    % relative uncertainty on the density
        neutron_uncert = 0.08;	% relative uncertainty on the local neutron production rate (mainly from scaling)
        muon_uncert = 0.08;		% relative uncertainty on the local production rate (mainly from scaling)
        % PRODUCTION RATES ---------------------------------------------- %
        P = DepthProd_Heisinger_Balco2017(0,lat,lon,elev,10,'calc'); % neglects tchanges in production right due to elevation changes during deposition
        neutron = P(1);
        muon_s = P(2);
        muon_f = P(3);

    case '36Cl'
        [data,~,~]  = xlsread('datarock_ROtest.xls',2); % load chemistry data
        [geodata,~,raw] = xlsread('datarock_ROtest.xls',3); % load geogrpahical data (lat,elevation)
        ind = input('Enter the number of the sample you want to run (the row within the data file) ');
        inds = ind;
        % check if there are more samples from this section in the data
        % file
        if length(find(geodata(:,7) == geodata(ind,7))) > 1
            inds = find(geodata(:,7) == geodata(ind,7));  % add indices of other files
        end

        % update constants
        % CONSTANTS ----------------------------------------------------- %
        rho = 2.2;              % density g/cmï¿½
        lambda = log(2)/3.013e5;% decay constant for 36Cl (Audi, 2017)
        PsCa = 52.16;           % Ca spallation rate, (Marrero et al., 2015)
        lat = geodata(ind,1);   % sample latitude
        elev = geodata(ind,2);  % sample elvation
        name = raw(ind+1,1);    % sample name
        age = geodata(inds,5)*1e3;    % sample age in yrs
        age_uncert = geodata(inds,6)*1e3;
        depth = geodata(inds,3)*100;  % samle depth
        sample_chem = data(ind+3,1:66);% sample chemistry
        % UNCERTAINTIES ------------------------------------------------- %
        rho_uncert = 0.1;	    % relative uncertainty on the density
        neutron_uncert = 0.08;	% relative uncertainty on the local neutron production rate (mainly from scaling)
        muon_uncert = 0.08;		% relative uncertainty on the local production rate (mainly from scaling)

        % assign scaling factor for nucleonic and muonic production ----- %
        pressure = stdatm(elev);        % convert to air pressure
        EL_f = stone2000(lat,pressure); % scaling neutrons
        EL_mu = Sc_muons(pressure);     % scaling muons
        % Calculate depth scalings
        so_e  = exp(-(1:depth(end)) .* rho/att_neut); % depth scaling neutrons
        so_mu = exp(-(1:depth(end)) .* rho/att_sm);   % depth scaling muons

        % Depth-scaled PRODUCTION RATES --------------------------------- %
        P_cosmo = nan(depth(end),1);  % this requires input to be sorted in a way that the younger samples are towards the top of the input file
        P_rad = nan(depth(end),1);
        % ignoring fast muons?!
        for i = 1:depth
            [P_cosmo(i),P_rad(i)] = clrock_mod(sample_chem,e,att_neut,so_e(i),so_mu(i),EL_f,EL_mu,PsCa,rho); % P_cosmo - cosmogenic and P_rad - radiogenic,
            % the above does not include attenuation length uncertainty
            % into the Production rate estimate. However, having this
            % function within the sam
        end
end

% Add additonal constraints on deposition ------------------------------- %
Model = struct();
Model.depths = [0];          % depths of constrained layers, first should always be 0 for surface (cm)
Model.depths_uncerts = [0];
Model.ages = [10e3];         % ages of constrained layers in yrs
Model.age_uncerts = [2e3];   % age uncertainty of constrained layers in yrs
% add constraints from main sample to additional ones
Model.depths = [Model.depths;depth];   depth = Model.depths;
Model.depths_uncerts = [Model.depths_uncerts; depth_uncert]; depth_uncert = Model.depths_uncerts;
Model.ages   = [Model.ages;age];       age   = Model.ages;
Model.age_uncerts = [Model.age_uncerts;age_uncert]; age_uncert = Model.age_uncerts;


%% START FORWARD MODELS ------------------------------------------------- %
% Final Storage matrix
Production = zeros(1,n);
% Initialize a progress bar to monitor the progression of the number of simulations n
wb = waitbar(0,'Welcome to the jungle...');

for i = 1:n
        % RANDOM SAMPLING ----------------------------------------------- %
        % Random sampling of the different parameters from a normal distribution for each of the n realisations of the simulation.
        % The ages of the core are sampled randomely within the uncertainty bounds of the dating provided to the function in Depth_age
        section_depth = round(normrnd(depth,depth_uncert));  % depth in cm
        Depth_age_guess = [section_depth,round(normrnd(age,age_uncert))];	% depth, age, matrix

        % if there's an age inversion, take a new sample
        counter = 0;
        while any(diff(Depth_age_guess(:,2)) < 0)
            section_depth = normrnd(depth,depth_uncert);  % depth in cm
            Depth_age_guess = [section_depth,round(normrnd(age,age_uncert))];
            counter = counter + 1;
            if counter > 1e4
                error("Couldn't sample a sequence without age inversion. Check yourage priors. ")
            end
        end

        % A column of densities is produced taking into account the uncertainty
        density_guess = normrnd(ones(1,section_depth)*rho, ones(1,section_depth)*rho*rho_uncert)';

        %For attenuation length we sample from a truncated (2 sigma on both sides) normal distribution to avoid negative and extreme values
        att_neut_guess = att_neut + att_neut* TruncatedGaussian(att_neut_uncert, [-2*att_neut_uncert, 2*att_neut_uncert],1);
        att_sm_guess = att_sm + att_sm* TruncatedGaussian(att_sm_uncert, [-2*att_sm_uncert, 2*att_sm_uncert],1);
        att_fm_guess = att_fm + att_fm* TruncatedGaussian(att_fm_uncert, [-2*att_fm_uncert, 2*att_fm_uncert],1);

        % assemble random production rates
        switch nuclide
            case '10Be'
                neutron_guess = normrnd(neutron, neutron*neutron_uncert);
                muon_s_guess = normrnd(muon_s, muon_s* muon_uncert);
                muon_f_guess = normrnd(muon_f, muon_f* muon_uncert);
            case '36Cl'
                PsCa_guess =  normrnd(PsCa, PsCa* neutron_uncert);   % uncertainty in Ca-spallation rate
                for ii = 1:depth(end)
                    avg_density = sum(density_guess(1:ii))/ii;
                    [P_cosmo(ii),P_rad(ii)] = clrock_mod(sample_chem,e,att_neut_guess,so_e(ii),so_mu(ii),EL_f,EL_mu,PsCa,avg_density); % P_cosmo - cosmogenic and P_rad - radiogenic
                end
        end
        %
        % First part of the routine
        % A random sediment layer sequence is built for each section of the core that is bound by two dates (minimum age and maximum age constrain)
        % This random sequence is stored in the matrix Sed_col, which contains layer thickness (Sed_col[,1]) and waiting time for the next layer = exposure time (Sed_col[,2])
        % Layer thickness and waiting times are drawn from two power law distributions as observed in many stockastic environments
        % For each section of the core the total length and time constrains are met (sum of length of each layer should be equal to total length of the core section and same for time)
        % In other words this means that for each section the average accumulation rate is the same as calculated from the core data
        % Each dated section of the core is filled with power law distributed layer thichnesses until the layer is full (dz = 0)
        Intervals = diff(Depth_age_guess);
        [r,~] = size(Intervals);
        Sed_col = [];
        for j  = 1:r

            counter = 0;
            tmp_thickness = [];
            while length(tmp_thickness) <= 2    % if a stratigraphy consists of only 2 or less layers, discard it and try again.

                tmp_thickness = [];
                dz = Intervals(j,1);
                while dz > 0            % This assembles a random layer stratigraphy
                     tmp = round(exp(exprnd(1/p_thickness)));   % this is a bit different in R, there the argument is rate = 1/mean
                     if (tmp > dz)
                         tmp = dz;
                     end
                     dz = dz - tmp;
                     tmp_thickness = [tmp_thickness;tmp];
                end
                counter = counter +1;
                if counter > 1e5
                    error(' stuck in the while loop for 1e5 iterations. Somethings wrong with your depth set up')
                end
            end

            % Now that the number of layers for the section is known, an equal number of waiting times is drawn also from a power law distribution
            % To make sure that the total amount of time is equal to the time covered by the section (i.e. the difference between lower and upper age) the waiting times are rescaled to that duration
            % (!! need to check if this is allowed, i.e. is the distribution of the rescalled waiting times also a power law distribution?)
            dt = Intervals(j,2);
            pos_time = exp(exprnd(1/p_time,length(tmp_thickness),1));
            tmp_time = dt*pos_time/sum(pos_time);

            % The data of this section is attached on top of the data from the previous section
            Sed_col = [[tmp_thickness,tmp_time] ; Sed_col]; % bind all data together
        end
%         Sed_col = Sed_col(1:end-1,:);       % why remove the last????

        % Second part of the routine
        % The postdepositional nuclide production can be computed for the simulated core based on the depth/age constrain imposed by the Sed_col matrix
        % The production of nuclides is computed from bottom to top (the loop indices decreases).
        % Each part of the loop calculates the nuclide production of a given layer as well as all layers under it.
        % Each iteration is then summed with the previous one
        profile = zeros(1,section_depth);
        for k =  length(Sed_col):-1:1
            % For each layer temporary depth, density, time and concentrtaion vectors are created
            % These vectors are made such that each element of the vector corresponds to 1 cm.
            tmp_depth = 1:sum(Sed_col(length(Sed_col):-1:k,1));
            tmp_density = cumsum(density_guess((section_depth-length(tmp_depth)+1):section_depth));
            tmp_time = repmat(Sed_col(k,2), length(tmp_depth),1);
            switch nuclide
                case '10Be'
                    tmp_conc = ((neutron_guess.*exp(-tmp_depth'./att_neut_guess))./lambda).*(1-exp(-tmp_time.*lambda))+...
                                ((muon_s_guess.*exp(-tmp_depth'./att_sm_guess))./lambda).*(1-exp(-tmp_time.*lambda))+...
                                ((muon_f_guess.*exp(-tmp_depth'./att_fm_guess))./lambda).*(1-exp(-tmp_time.*lambda));
                case '36Cl'
                    tmp_conc = ((P_cosmo(tmp_depth).*exp(-tmp_depth'./att_neut_guess))./lambda).*(1-exp(-tmp_time.*lambda))+ P_rad(1);
            end

            % The result is added to the production vector
            profile((section_depth-length(tmp_depth)+1):section_depth) = profile((section_depth-length(tmp_depth)+1):section_depth) + tmp_conc;
        end
        Production(i) = profile(end);
    % Increment progress bar
    waitbar(i/(n),wb)
end
close(wb)

Pmean = mean(Production);
Pmedian = median(Production);
Pstd = std(Production);
% Pquantiles = quantile(Production,[0.022,0.977],2);
disp(['postburial production = ' num2str(Pmean) ' +/- ' num2str(Pstd)]) 

%% EXPORT THE DATA
if export
    filename = input('What name do you want to give this run for saving? ');
%     print(['PostburialProd_' name '_' filename '.eps'], '-depsc', '-painters') % save figure
    savefig(h,['./output/PostburialProd_' ,name{1} ,'_' ,filename, '.fig'])

    % add all parameters to model
    Model.density = rho;
    Model.att_neut = att_neut;
    Model.att_fm  = att_fm;
    Model.att_sm = att_sm;
    Model.p_thickness = p_thickness;
    Model.p_time = p_time;
    Model.nruns = n;
    Model.mean_postburial = Pmean;
    Model.sd_postburial = Pstd;

    save(['./output/PostburialProd_' name{1} '_' num2str(filename) '.mat'], 'model')                % save model parameters
end
