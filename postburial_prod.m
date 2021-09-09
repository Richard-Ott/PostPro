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

%%%%%%%%%%%%%%%%%%%%%%%
% I think the power law formulation is incorrect. Gomez uses a log base of
% 2 and not the ln. Also, this cide takes the exponential of exprnd and I
% think thats incorrect because you take the exp of the already exponenated
% number
%%%%%%%%%%%%%%%%%%%%%%

% Richard Ott, 2020
clc
clear
close all

addpath '.\subroutines'
addpath '..\..\..\..\Crete\Cretan_fans\data'

% USER CHOICE ----------------------------------------------------------- %
nuclide = '10Be';       % Choose '10Be' or '36Cl'
export = 0;             % do you want to save figures and model data
n = 1e3;                % number of runs
global scaling_model
scaling_model = 'lm';   % choose your scaling model, nomenclature follows Cronus

% load sample data in Cronus excel format
% [num,txt,~] = xlsread('36Cl_data_CRONUS.xlsx','Matlab Postburial');
[num,txt,~] = xlsread('10Be_data_CRONUS','Matlab Postburial');

%% assign data and constants -------------------------------------------- %
p_thickness = 2.06; 	% parameter for power law distribution of sediment layer thickness (Gomez et al., 2002) 2.06
p_time = 1.4;		    % parameter for power law distribution of sediment layer waiting time (Gomez et al., 2002)
e = 0;                  % erosion rate (mm/yr)

% Production rate uncertainties
switch nuclide
    case '10Be'
        Ps_uncert = 0.08;       % uncertainty spallation production 10Be (Phillips et al., 2015)
        Pmu_uncert = 0.3;       % uncertainty muon production 10Be (Phillips 2015 has 10% uncert on f*, I triple this)
    case '36Cl'
        Ps_uncert  = 0.065;	    % relative uncertainty on the spallation production rate
        % average of PsCa and PsK (Cronus v2.1)
        Pmu_uncert = 0.25;      % relative uncertainty on the muon production rate
        % I am using the middle between the f*Ca and f*K uncertainty in Marrero, 2016
        Pth_uncert  = 0.36; % relative uncertainty on the thermal neutron capture production rate (Cronus v2.1 with lm)
        Peth_uncert = 0.35; % relative uncertainty on the thermal neutron capture production rate (Cronus v2.1 with lm)
end

% Add additonal constraints on deposition ------------------------------- %
% here you need to enter the parameter of the layers that are not
% constrained by the sample you are running. This includes the surface and
% its age as well as additional layers in you section that may be dated and
% therefore constrain the deposition of overburden for your sample
Model = struct();
Model.depths         = [0];  % depths of constrained layers (cm), first should always be 0 for surface (cm)
Model.depths_uncerts = [0];  % (cm)
Model.ages           = [71.1e3];    % ages of constrained layers in yrs
Model.age_uncerts    = [7.5e3];  % age uncertainty of constrained layers in yrs
% Model.comment        = 'The age of the top of deposit is taken as the middle between the age of the higher OSL sample and the C14 date from river incision the uncertainty corresponds to the difference between the middle value and the OSL/14C age'; 

%% Load data                       
switch nuclide
    case '10Be'
        ind  = input('Enter the number of the sample you want to run (the row within the data file) ');
        depth= num(ind,14)*100;        % depth of sample (cm), I enter the values as m into the excel sheet even though the Cronus input is in g/cmÂ²
        depth_uncert = num(ind,29)*100;% uncertainty in overburden (cm)
        elev = num(ind,3);             % elevation in (m)
        Tshd = num(ind,7);             % topographic shielding factor
        name = txt(ind,1);             % sample name
        age = num(ind,33)*1e3;         % age of sample   (yrs)
        age_uncert = num(ind,34)*1e3;  % age uncertainty (yrs)
        rho = num(ind,6);              % density of burying material n g/cm³

        % CONSTANTS ----------------------------------------------------- %
        lambda = log(2)/1.39e6;	       % decay constant for 10Be
        num(ind,4) = stdatm(elev);     % convert to air pressure
        
    case '36Cl'

        ind = input('Enter the number of the sample you want to run (the row within the data file) ');

        % CONSTANTS ----------------------------------------------------- %
        lambda = log(2)/3.013e5;       % decay constant for 36Cl (Audi, 2017)
        elev = num(ind,3);             % sample elvation
        name = txt(ind,1);             % sample name
        age = num(ind,80)*1e3;         % sample age in yrs
        age_uncert = num(ind,81)*1e3;
        depth = num(ind,12)*100;       % sample depth in cm
        depth_uncert = num(ind,51)*100;
        rho = num(hallo);              % density of burying material in g/cm³

        % assign scaling factor for nucleonic and muonic production ----- %
        num(ind,4) = stdatm(elev);    % convert to air pressure
end

% add constraints from main sample to additional ones
Model.depths         = [Model.depths;depth];  depth = Model.depths; 
Model.depths_uncerts = [Model.depths_uncerts; depth_uncert]; depth_uncert = Model.depths_uncerts; 
Model.ages           = [Model.ages;age];      age   = Model.ages;
Model.age_uncerts    = [Model.age_uncerts;age_uncert];       age_uncert = Model.age_uncerts;


%% PRODUCTION RATES ----------------------------------------------------- %
pp=physpars();
max_depth = depth(end)*rho + 5*depth_uncert(end)*rho; % get maximum possible depth to come up in forward models, this will be the depth until which the production profile is calculated, in g/cm²
max_age   = age(end)   + 5*age_uncert(end);   % get maximum possible age for forward models, 4-sigma cutoff
num = num(ind,:);  % I think this can be removed?!?!
switch nuclide
    case '10Be'
        num([14,29]) = 0;                     % set "depth-to-top-of-sample" for CRONUS back to zero in order to calculate full production rate profile
        [nominal10,uncerts10] = createage1026(num);               % get basic sample info 
        sp = samppars1026(nominal10);                             % Extract the sample parameters from the sampledatavector.
        sf = scalefacs1026(sp,scaling_model);                     % get scaling factors
        cp = comppars1026(pp,sp,sf,max_depth);                    % computed parameters
        
        % for all potential depths and ages calculate the production rate
        % profile. 
        % Production rates are sorted in matrix with depth rows (in g/cm2) and time
        % columns (every 100 years)
        Ps10  = nan(max_depth,ceil(max_age/1e2));
        Pmu10 = nan(max_depth,ceil(max_age/1e2));
        for i = 1:size(Ps10,2)
            % this scales back the productoin rate to a certain time, age has to be negative and in (ka)
            sf.currentsf=getcurrentsf(sf,-i*1e-1,scaling_model,'be');  
            % get production rates for a production profile
            [~,~,Ps10(:,i),Pmu10(:,i),~,~] = prodz1026(1:max_depth,pp,sf,cp);   
        end
        
    case '36Cl'
        num([12,51]) = 0;      % set "depth-to-top-of-sample" for CRONUS back to zero in order to calculate full production rate profile
        [nominal36,uncerts36,cov36]=createage36(num);    % Get basic info about the sample ages.
        sp = samppars36(nominal36);         % Extract the sample parameters from the sampledatavector.
        sf = scalefacs36(sp,scaling_model);  % get scaling factors
        cp = comppars36(pp,sp,sf,max_depth); % computed parameters
        
        % for all potential depths and ages calculate the production rate
        % profile. 
        % Production rates are sorted in matrix with depth rows (in g/cm2) and time
        % columns (every 100 years)
        Ps36   = nan(max_depth,ceil(max_age/1e2));
        Pmu36  = nan(max_depth,ceil(max_age/1e2));
        Pth36  = nan(max_depth,ceil(max_age/1e2));
        Peth36 = nan(max_depth,ceil(max_age/1e2));
        for i = 1:size(Ps36,2)
            % this scales back the production rate to a certain time, age has to be negative and in (ka)
            sf.currentsf=getcurrentsf(sf,-i*1e-1,scaling_model,'cl');     
            % get production rates for a production profile
            [~,Ps36(:,i),~,~,~,~,Pth36(:,i),Peth36(:,i),Pmu36(:,i),~,~,~,~,~,~,~] = ...
            prodz36(1:max_depth,pp,sf,cp);             
        end
end


%% START FORWARD MODELS ------------------------------------------------- %
% Final Storage matrix for postburial production at sample depth
PostProduction = zeros(1,n);

% draw random production rates before loop to speed up computation,
% truncate normal distribution to avoid negative production rates
Ps_rand  = truncnormrnd([n,1],0,Ps_uncert ,-1,1);
Pmu_rand = truncnormrnd([n,1],0,Pmu_uncert,-1,1);
if strcmpi(nuclide,'36Cl')
    Pth_rand  = truncnormrnd([n,1],0,Pth_uncert ,-1,1);
    Peth_rand = truncnormrnd([n,1],0,Peth_uncert,-1,1);
end
        

% Initialize a progress bar to monitor the progression of the number of simulations n
wb = waitbar(0,'Welcome to the jungle...');

for i = 1:n
        % RANDOM SAMPLING ----------------------------------------------- %
        % Random sampling of the different parameters from a normal distribution 
        % for each of the n realisations of the simulation. The ages of the
        % core are sampled randomely within the uncertainty bounds of the 
        % dating provided to the function in Depth_age to avoid computational
        % problems the normal distribution gets truncated at 3 sigma
        section_depth = round(truncnormrnd(1,depth(end),depth_uncert(end),depth(end)-3*depth_uncert(end), depth(end)+3*depth_uncert(end)));  % depth in cm 
        Depth_age_guess = [[depth(1:end-1);section_depth],round(normrnd(age,age_uncert))];	% depth, age, matrix

        % if there's an age inversion, take a new sample
        counter = 0;
        while any(diff(Depth_age_guess(:,2)) <= 1) || any(Depth_age_guess(:,2)<0) || any(diff(Depth_age_guess(:,1)) <= 5)
            section_depth = round(truncnormrnd(1,depth(end),depth_uncert(end),depth(end)-3*depth_uncert(end), depth(end)+3*depth_uncert(end)));  % depth in cm
            Depth_age_guess = [[depth(1:end-1);section_depth],round(normrnd(age,age_uncert))];
            counter = counter + 1;
            if counter > 1e4
                error("Couldn't sample a sequence without age inversion. Check your age priors. ")
            end
        end

        % assemble random production rates
        switch nuclide
            case '10Be'
                Ps  = Ps10  + Ps10  .* Ps_rand(i);
                Pmu = Pmu10 + Pmu10 .* Pmu_rand(i);
            case '36Cl'
                Ps   = Ps36  + Ps36    .* Ps_rand(i);
                Pmu  = Pmu36 + Pmu36   .* Pmu_rand(i);
                Pth  = Pth36  + Pth36  .* Pth_rand(i);
                Peth = Peth36 + Peth36 .* Peth_rand(i);               
        end
        %
        % First part of the routine
        % A random sediment layer sequence is built for each section of the 
        % core that is bound by two dates (minimum age and maximum age 
        % constrain). This random sequence is stored in the matrix Sed_col,
        % which contains layer thickness (Sed_col[,1]) and waiting time for
        % the next layer = exposure time (Sed_col[,2]). Layer thickness and
        % waiting times are drawn from two power law distributions as 
        % observed in many stockastic environments. For each section of the 
        % core the total length and time constrains are met (sum of length
        % of each layer should be equal to total length of the core section 
        % and same for time). In other words this means that for each 
        % section the average accumulation rate is the same as calculated 
        % from the core data. Each dated section of the core is filled with
        % power law distributed layer thichnesses until the layer is full (dz = 0)
        
        Intervals = diff(Depth_age_guess);
        [r,~] = size(Intervals);
        Sed_col = [];
        for j  = 1:r
            newtry = 0;
            while newtry == 0 % this while loop ensures that if a weird sample of thicknesses is drawn the 
                % algorithm does not get stuck trying to find a suitable
                % age distribution
            
                counter = 0;
                tmp_thickness = [];
                while length(tmp_thickness) <= 2    % if a stratigraphy consists of only 2 or less layers, discard it and try again.

                    tmp_thickness = [];
                    dz = Intervals(j,1);
                    while dz > 0            % This assembles a random layer stratigraphy
                         tmp = round(exp(exprnd(1/p_thickness)));  
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

                % Now that the number of layers for the section is known, an 
                % equal number of waiting times is drawn also from a power law 
                % distribution. To make sure that the total amount of time is
                % equal to the time covered by the section (i.e. the difference
                % between lower and upper age) the waiting times are being
                % drawn until the the correct total duration (+-10% is
                % achieved). This is necessary because there is no exponential
                % sampling function with fixed sum for matlab. The waiting times
                % with +-20% difference are then rescaled to the exact total time.
                dt = Intervals(j,2);
                pos_time = exp(exprnd(1/p_time,length(tmp_thickness),1)); % draw random times
                counter = 0;
                newtry = 1;
                while abs(sum(pos_time) - dt) > dt*0.2  && newtry ==0 % difference sampled time to total time larger 10%, then resample
                    pos_time = exp(exprnd(1/p_time,length(tmp_thickness),1)); % draw random times
                    counter = counter+1;
                    if counter > 1e3                 % if the algorithm cannot find a suiting age distribution draw a new thcikness distribution
                        newtry = 0;
                    end
                end
                tmp_time = dt*pos_time/sum(pos_time);   % rescale to total time dt

                % The data of this section is attached on top of the data from the previous section
                Sed_col = [[tmp_thickness,tmp_time] ; Sed_col]; % bind all data together
            end
        end

        % Second part of the routine
        % The postdepositional nuclide production can be computed for the 
        % simulated core based on the depth/age constrain imposed by the 
        % Sed_col matrix.
        % The production of nuclides is computed from bottom to top (the 
        % loop indices decreases). Each part of the loop calculates the 
        % nuclide production of a given layer as well as all layers under it.
        % Each iteration is then summed with the previous one
        Sed_col_gcm2 = Sed_col;  Sed_col_gcm2(:,1) = round(Sed_col(:,1).*rho);  % convert depth to gm/2 to match production rates
        section_depth_gcm2 = round(section_depth * rho);  % convert to g/cm2 to match production rate tables
        profile = zeros(1,section_depth_gcm2);
        for k =  length(Sed_col):-1:1
            % For each layer temporary depth, density, time and concentration vectors are created
            % These vectors are made such that each element of the vector corresponds to 1 cm.
            tmp_depth = 1:sum(Sed_col_gcm2(length(Sed_col_gcm2):-1:k,1));
            tmp_time  = repmat(Sed_col_gcm2(k,2), length(tmp_depth),1);
            
            % get mean depositional age of this layer for correction
            % scaling factor
            if k == length(Sed_col_gcm2)
                tmp_age   = Depth_age_guess(end,2) - cumsum(Sed_col_gcm2(k:length(Sed_col_gcm2),2))/2;
            else
                tmp_age   = Depth_age_guess(end,2) - (sum(Sed_col_gcm2(k:length(Sed_col_gcm2),2))+sum(Sed_col_gcm2(k+1:length(Sed_col_gcm2),2)))/2;
            end
            
            Pind = round(tmp_age/1e2)+1;  % get index of production rate column that is closest in age to tmp_age
            if Pind > round(max_age/1e2)
                error('The random sample is older than the oldest computed production rate. Increase the safety factor for max_age or use a truncated normal distribution for drawing the samples')
            end
            switch nuclide
                case '10Be'
                    tmp_conc = (Ps(tmp_depth,Pind) + Pmu(tmp_depth,Pind)).*tmp_time;  % Production at/g
                    tmp_conc = tmp_conc.* exp(-tmp_time.*lambda);      % radioactive decay
                case '36Cl'
                    tmp_conc = (Ps(tmp_depth,Pind) + Pmu(tmp_depth,Pind) + ...
                        Peth(tmp_depth,Pind) + Pth(tmp_depth,Pind)).*tmp_time;       % Production at/g
                    tmp_conc = tmp_conc.* exp(-tmp_time.*lambda);      % radioactive decay
            end

            % The result is added to the production vector
            profile((section_depth_gcm2-length(tmp_depth)+1):section_depth_gcm2) = profile((section_depth_gcm2-length(tmp_depth)+1):section_depth_gcm2) + tmp_conc';
        end
        
        % add production after end of deposition
        if Depth_age_guess(1,2) ~= 0 % if there is a no-deposition period at the end
            Pind = round(Depth_age_guess(1,2)/1e2);
            switch nuclide
                case '10Be'
                    profile = profile + (sum(Ps(1:section_depth_gcm2,1:Pind) + ...
                        Pmu(1:section_depth_gcm2,1:Pind),2).*1e2);        % Production at/g
                case '36Cl'
                    profile = profile + (sum(Ps(1:section_depth_gcm2,1:Pind) + ...
                        Pmu(1:section_depth_gcm2,1:Pind) + Pth(1:section_depth_gcm2,1:Pind) ...
                    + Peth(1:section_depth_gcm2,1:Pind),2).*1e2);        % Production at/g
            end
            profile = profile.* exp(-1e2.*lambda);    % radioactive decay
        end
        
        PostProduction(i) = profile(end);
    % Increment progress bar
    waitbar(i/(n),wb)
end
close(wb)

Pmean = round(mean(PostProduction));
Pmedian = round(median(PostProduction));
Pstd = round(std(PostProduction));
Pquartiles = round(quantile(PostProduction,[0.25,0.75]));
Pquantiles = round(quantile(PostProduction,[0.17,0.83]));
disp(['postburial production ' name{1} ' = ' num2str(Pmean) ' +/- ' num2str(Pstd)]) 
disp(['postburial production median ' name{1} ' = ' num2str(Pmedian) ' +/- ' num2str(diff(Pquantiles)/2)]) 
disp(['postburial production quantiles ' name{1} ' = ' num2str(Pquantiles(1)) ' - ' num2str(Pquantiles(2))]) 


%% EXPORT THE DATA
if export
    filename = input('What name do you want to give this run for saving? ','s');

    % add all parameters to model
    Model.scaling_model = scaling_model;
    Model.Ps_uncert   = Ps_uncert;
    Model.Pmu_uncert  = Pmu_uncert;
    Model.p_thickness = p_thickness;
    Model.p_time = p_time;
    Model.nruns = n;
    Model.mean_postburial = Pmean;
    Model.sd_postburial = Pstd;
    Model.median = Pmedian;
    Model.quartiles = Pquartiles;
    Model.quants1783 = Pquantiles;

    save(['./output/' nuclide '/PostburialProd_' name{1} '_' num2str(filename) '.mat'], 'Model')                % save model parameters
end
