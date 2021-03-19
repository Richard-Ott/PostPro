% This code calculates depth profiles of post burial production for 10Be
% and 36Cl (could easily expanded to work for 26Al and 14C). The
% calculations are based on geochronologic anchor ages. For the sites with
% geochronologic ages the site specific parameters (latitude, elevation,
% shielding, if applicable chemistry etc) need to be specified. The code
% then takes the one or more samples and runs Monte Carlo simulation of
% deposition, where deposition of a random amount of layers with certain
% thickness and time intervals follow power law distributions, as suggested
% by Gomez et al. 2002.
% The code in this branch uses CRONUS to determine production rates.
% The code is based on an R-script by M. Lupker but is now heavily
% modified.
%
% Richard Ott, 2020
clc
clear
close all

addpath '.\subroutines'
% addpath '..\data'
addpath '..\..\..\..\Crete\Cretan_fans\data'

% USER CHOICE ----------------------------------------------------------- %
nuclide = '10Be';       % Choose '10Be' or '36Cl'
export = 1;             % do you want to save figures and model data
n = 5e3;                % number of runs
scaling_model = 'st';   % choose your scaling model, nomenclature follows Cronus

% load sample data in Cronus excel format
[num,txt,~] = xlsread('Cl_Crete_617_forCronusCalc.xlsx','sheet2');

%% assign data and constants -------------------------------------------- %
p_thickness = 2.06; 	% parameter for power law distribution of sediment layer thickness (Gomez et al., 2002)
p_time = 1.4;		    % parameter for power law distribution of sediment layer waiting time (Gomez et al., 2002)
e = 0;                  % erosion rate (mm/yr)
P_neut_uncert = 0.05;	% relative uncertainty on the spallation production rate
P_sm_uncert   = 0.3;    % relative uncertainty on the stopping muons production rate
P_fm_uncert   = 0.3;    % relative uncertainty on the fast muons production rate
if strcmpi(nuclide,'36Cl')
    P_th_uncert  = 0.5; % relative uncertainty on the thermal neutron capture production rate
    P_eth_uncert = 0.5; % relative uncertainty on the thermal neutron capture production rate
end

% Add additonal constraints on deposition ------------------------------- %
Model = struct();
Model.depths = [0; 500; 550];          % depths of constrained layers, first should always be 0 for surface (cm)
Model.depths_uncerts = [0; 100; 100];
Model.ages = [8e3; 28e3; 49e3];         % ages of constrained layers in yrs
Model.age_uncerts = [2.5e3; 2e3; 4e3];   % age uncertainty of constrained layers in yrs
% add constraints from main sample to additional ones
Model.depths = [Model.depths;depth];   depth = Model.depths;
Model.depths_uncerts = [Model.depths_uncerts; depth_uncert]; depth_uncert = Model.depths_uncerts;
Model.ages   = [Model.ages;age];       age   = Model.ages;
Model.age_uncerts = [Model.age_uncerts;age_uncert]; age_uncert = Model.age_uncerts;

switch nuclide
    case '10Be'
        ind = input('Enter the number of the sample you want to run (the row within the data file) ');
        depth = data(ind,11)*100;        % depth of sample (cm)
        depth_uncert = data(ind,12)*100; % uncertainty in overburden (cm)
        age = data(ind,13)*1e3;         % age of sample   (yrs)
        age_uncert = data(ind,14)*1e3;  % age uncertainty (yrs)
        lat = data(ind,3);              % latitude in degree, must be same for all depths
        lon = data(ind,4);              % longitude in degree (needed for cutoff rigidity)
        elev = data(ind,5);             % elevation in (m)
        Tshd = data(ind,15);            % topographic shielding factor
        name = raw(1+ind,1);

        % CONSTANTS ----------------------------------------------------- %
        lambda = log(2)/1.39e6;	% decay constant for 10Be
        
    case '36Cl'

        ind = input('Enter the number of the sample you want to run (the row within the data file) ');
        inds = ind;
        % check if there are more samples from this section in the data
        % file
        if length(find(geodata(:,7) == geodata(ind,7))) > 1
            inds = find(geodata(:,7) == geodata(ind,7));  % add indices of other files
        end

        % update constants
        % CONSTANTS ----------------------------------------------------- %
        lambda = log(2)/3.013e5;% decay constant for 36Cl (Audi, 2017)
        lat = geodata(ind,1);   % sample latitude
        elev = geodata(ind,2);  % sample elvation
        name = raw(ind+1,1);    % sample name
        age = geodata(inds,5)*1e3;    % sample age in yrs
        age_uncert = geodata(inds,6)*1e3;
        depth = geodata(inds,3)*100;  % samle depth
        sample_chem = data(ind+3,1:66);% sample chemistry


        % assign scaling factor for nucleonic and muonic production ----- %
        pressure = stdatm(elev);        % convert to air pressure
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

%% PRODUCTION RATES ----------------------------------------------------- %
pp=physpars();
max_depth = depth(end) + 10*depth_uncert(end); % get maximum possible depth to come up in forward models, this will e the depth until which the production profile is calculated
switch nuclide
    case '10Be'
        sp = samppars1026(sampledata);        % Extract the sample parameters from the sampledatavector.
        sf = scalefacs1026(sp,scaling_model); % get scaling factors
        cp = comppars1026(pp,sp,sf,maxdepth); % computed parameters
        [~,~,Ps10,Pmu10,~,~] = prodz1026(1:max_depth,pp,sf,cp);   % get production rates for a productoin profile
    case '36Cl'
        sp = samppars36(sampledata);        % Extract the sample parameters from the sampledatavector.
        sf = scalefacs36(sp,scaling_model); % get scaling factors
        cp = comppars36(pp,sp,sf,maxdepth); % computed parameters
        [~,Ps36,~,~,~,~,Pth36,Peth36,Pmu36,~,~,~,~,~,~,~] = ...
        prodz36(currentdepths-deltadepth/2,pp,sf,cp);             % get production rates for a productoion profile
end


%% START FORWARD MODELS ------------------------------------------------- %
% Final Storage matrix
Production = zeros(1,n);
% Initialize a progress bar to monitor the progression of the number of simulations n
wb = waitbar(0,'Welcome to the jungle...');

for i = 1:n
        % RANDOM SAMPLING ----------------------------------------------- %
        % Random sampling of the different parameters from a normal distribution for each of the n realisations of the simulation.
        % The ages of the core are sampled randomely within the uncertainty bounds of the dating provided to the function in Depth_age
        section_depth = round(truncnormrnd(1,depth(end),depth_uncert(end),600 , 1000));  % depth in cm
%         section_depth = round(normrnd(depth(end),depth_uncert(end)));  % depth in cm
        Depth_age_guess = [[depth(1:end-1);section_depth],round(normrnd(age,age_uncert))];	% depth, age, matrix

        % if there's an age inversion, take a new sample
        counter = 0;
        while any(diff(Depth_age_guess(:,2)) < 0)
            section_depth = round(truncnormrnd(1,depth(end),depth_uncert(end),600 , 1000));  % depth in cm
%             section_depth = round(normrnd(depth(end),depth_uncert(end)));  % depth in cm
            Depth_age_guess = [[depth(1:end-1);section_depth],round(normrnd(age,age_uncert))];
            counter = counter + 1;
            if counter > 1e4
                error("Couldn't sample a sequence without age inversion. Check yourage priors. ")
            end
        end


        %For production rates we sample from a truncated (3 sigma on both sides) normal distribution to avoid negative and extreme values
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
                    tmp_conc = ((P_cosmo(tmp_depth).*exp(-tmp_depth'./att_neut_guess))./lambda).*(1-exp(-tmp_time.*lambda))+ P_rad(1); % the radiogenic production is incorrect. it also needs to be multiplied by time
            end

            % The result is added to the production vector
            profile((section_depth-length(tmp_depth)+1):section_depth) = profile((section_depth-length(tmp_depth)+1):section_depth) + tmp_conc';
        end
        Production(i) = profile(end);
    % Increment progress bar
    waitbar(i/(n),wb)
end
close(wb)

Pmean = round(mean(Production));
Pmedian = round(median(Production));
Pstd = round(std(Production));
Pquantiles = round(quantile(Production,[0.17,0.83]));
disp(['postburial production ' name{1} ' = ' num2str(Pmean) ' +/- ' num2str(Pstd)]) 
disp(['postburial production median ' name{1} ' = ' num2str(Pmedian) ' +/- ' num2str(diff(Pquantiles)/2)]) 
disp(['postburial production quantiles ' name{1} ' = ' num2str(Pquantiles(1)) ' - ' num2str(Pquantiles(2))]) 


%% EXPORT THE DATA
if export
    filename = input('What name do you want to give this run for saving? ','s');
%     print(['PostburialProd_' name '_' filename '.eps'], '-depsc', '-painters') % save figure
%     savefig(h,['./output/PostburialProd_' ,name{1} ,'_' ,filename, '.fig'])

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
    Model.median = Pmedian;
    Model.quants1783 = Pquantiles;

    save(['./output/PostburialProd_' name{1} '_' num2str(filename) '.mat'], 'Model')                % save model parameters
end
