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

% MAKE ASURE THAT UNITS OF G/CM² COMING OUT OF CRONUS GO WITH MY TREATING
% OF EVERYTHING  AS LENGTH UNITS...


% Richard Ott, 2020
clc
clear
close all

addpath '.\subroutines'
% addpath '..\data'
addpath '..\..\..\..\Crete\Cretan_fans\data'

% USER CHOICE ----------------------------------------------------------- %
nuclide = '10Be';       % Choose '10Be' or '36Cl'
export = 0;             % do you want to save figures and model data
n = 5e3;                % number of runs
scaling_model = 'sf';   % choose your scaling model, nomenclature follows Cronus

% load sample data in Cronus excel format
% [num,txt,~] = xlsread('36Cl_data_CRONUS.xlsx','sheet2');
[num,txt,~] = xlsread('10Be_data_CRONUS',2);

%% assign data and constants -------------------------------------------- %
p_thickness = 2.06; 	% parameter for power law distribution of sediment layer thickness (Gomez et al., 2002)
p_time = 1.4;		    % parameter for power law distribution of sediment layer waiting time (Gomez et al., 2002)
e = 0;                  % erosion rate (mm/yr)
Ps_uncert  = 0.05;	% relative uncertainty on the spallation production rate
Pmu_uncert = 0.3;    % relative uncertainty on the muon production rate
if strcmpi(nuclide,'36Cl')
    Pth_uncert  = 0.5; % relative uncertainty on the thermal neutron capture production rate
    Peth_uncert = 0.5; % relative uncertainty on the thermal neutron capture production rate
end

% Add additonal constraints on deposition ------------------------------- %
% here you need to enter the parameter of the layers that are not
% constrained by the sample you are running. This includes the surface and
% its age as well as additional layers in you section that may be dated and
% therefore constrain the deposition of overbrudn for your sample
Model = struct();
Model.depths         = [0];  % depths of constrained layers, first should always be 0 for surface (cm)
Model.depths_uncerts = [0];
Model.ages           = [55e3];  % ages of constrained layers in yrs
Model.age_uncerts    = [5e3];  % age uncertainty of constrained layers in yrs

switch nuclide
    case '10Be'
        ind  = input('Enter the number of the sample you want to run (the row within the data file) ');
        depth= num(ind,14)*100;        % depth of sample (cm), I enter the values as m into the excel sheet even though the Cronus input is in g/cmÂ²
        depth_uncert = num(ind,29)*100;% uncertainty in overburden (cm)
        lat  = num(ind,1);             % latitude in degree, must be same for all depths
        lon  = num(ind,2);             % longitude in degree (needed for cutoff rigidity)
        elev = num(ind,3);             % elevation in (m)
        Tshd = num(ind,7);            % topographic shielding factor
        name = txt(4+ind,1);
        
        % from separate input sheet
        age = num(ind,33)*1e3;         % age of sample   (yrs)
        age_uncert = num(ind,34)*1e3;  % age uncertainty (yrs)

        % CONSTANTS ----------------------------------------------------- %
        lambda = log(2)/1.39e6;	          % decay constant for 10Be
        num(ind,4) = stdatm(elev);        % convert to air pressure
        
    case '36Cl'

        ind = input('Enter the number of the sample you want to run (the row within the data file) ');
        inds = ind;

        % CONSTANTS ----------------------------------------------------- %
        lambda = log(2)/3.013e5;% decay constant for 36Cl (Audi, 2017)
        lat  = num(ind,1);      % sample latitude
        elev = num(ind,2);      % sample elvation
        name = txt(ind+1,1);    % sample name
        age = num(inds,5)*1e3;    % sample age in yrs
        age_uncert = num(inds,6)*1e3;
        depth = num(inds,3)*100;  % samle depth

        % assign scaling factor for nucleonic and muonic production ----- %
        num(XXXX) = stdatm(elev);        % convert to air pressure

end

% add constraints from main sample to additional ones
Model.depths = [Model.depths;depth];   depth = Model.depths;
Model.depths_uncerts = [Model.depths_uncerts; depth_uncert]; depth_uncert = Model.depths_uncerts;
Model.ages   = [Model.ages;age];       age   = Model.ages;
Model.age_uncerts = [Model.age_uncerts;age_uncert]; age_uncert = Model.age_uncerts;

%% PRODUCTION RATES ----------------------------------------------------- %
pp=physpars();
max_depth = depth(end) + 4*depth_uncert(end); % get maximum possible depth to come up in forward models, this will e the depth until which the production profile is calculated
max_age   = age(end)   + 4*age_uncert(end);   % get maximum possible age for forward models, 4-sigma cutoff
num = num(ind,1:end-2);
num([14,19]) = 0;
switch nuclide
    case '10Be'
        [nominal10,uncerts10] = createage1026(num);               % get basic sample info 
        sp = samppars1026(nominal10);                             % Extract the sample parameters from the sampledatavector.
        sf = scalefacs1026(sp,scaling_model);                     % get scaling factors
        cp = comppars1026(pp,sp,sf,max_depth);                     % computed parameters
        
        % for all potential depths and ages calculate the production rate
        % profile. 
        % Production rates are sorted in matrix with depth rows (in cm) and time
        % columns (every 100 years)
        Ps10  = nan(max_depth,ceil(max_age/1e2));
        Pmu10 = nan(max_depth,ceil(max_age/1e2));
        for i = 1:size(Ps10,2)
            sf.currentsf=getcurrentsf(sf,-i*1e-1,scaling_model,'be');     % this scales back the productoin rate to a certain time, age has to be negative and in (ka)
            [~,~,Ps10(:,i),Pmu10(:,i),~,~] = prodz1026(1:max_depth,pp,sf,cp);   % get production rates for a production profile
        end
        
    case '36Cl'
        [nominal36,uncerts36,cov36]=createage36(num);    % Get basic info about the sample ages.
        sp = samppars36(sampledata);        % Extract the sample parameters from the sampledatavector.
        sf = scalefacs36(sp,scaling_model); % get scaling factors
        cp = comppars36(pp,sp,sf,max_depth); % computed parameters
        
        % for all potential depths and ages calculate the production rate
        % profile. 
        % Production rates are sorted in matrix with depth rows (in cm) and time
        % columns (every 100 years)
        Ps36   = nan(max_depth,ceil(max_age/1e2));
        Pmu36  = nan(max_depth,ceil(max_age/1e2));
        Pth36  = nan(max_depth,ceil(max_age/1e2));
        Peth36 = nan(max_depth,ceil(max_age/1e2));
        for i = 1:size(Ps36,2)
            sf.currentsf=getcurrentsf(sf,-i*1e-1,scaling_model,'cl');     % this scales back the productoin rate to a certain time, age has to be negative and in (ka)
            [~,Ps36(:,i),~,~,~,~,Pth36(:,i),Peth36(:,i),Pmu36(:,i),~,~,~,~,~,~,~] = ...
            prodz36(1:max_depth,pp,sf,cp);             % get production rates for a production profile
        end
end


%% START FORWARD MODELS ------------------------------------------------- %
% Final Storage matrix
PostProduction = zeros(1,n);
% Initialize a progress bar to monitor the progression of the number of simulations n
wb = waitbar(0,'Welcome to the jungle...');

for i = 1:n
        % RANDOM SAMPLING ----------------------------------------------- %
        % Random sampling of the different parameters from a normal distribution for each of the n realisations of the simulation.
        % The ages of the core are sampled randomely within the uncertainty bounds of the dating provided to the function in Depth_age
        % to avoid computational problems the normal distribution gets
        % truncated at 4 sigma
        section_depth = round(truncnormrnd(1,depth(end),depth_uncert(end),depth(end)-4*depth_uncert(end), depth(end)+4*depth_uncert(end)));  % depth in cm
        Depth_age_guess = [[depth(1:end-1);section_depth],round(normrnd(age,age_uncert))];	% depth, age, matrix

        % if there's an age inversion, take a new sample
        counter = 0;
        while any(diff(Depth_age_guess(:,2)) < 0)
            section_depth = round(truncnormrnd(1,depth(end),depth_uncert(end),depth(end)-4*depth_uncert(end), depth(end)+4*depth_uncert(end)));  % depth in cm
            Depth_age_guess = [[depth(1:end-1);section_depth],round(normrnd(age,age_uncert))];
            counter = counter + 1;
            if counter > 1e4
                error("Couldn't sample a sequence without age inversion. Check yourage priors. ")
            end
        end

        % For production rates we sample from a truncated (3 sigma on both sides) normal distribution to avoid negative and extreme values
        % assemble random production rates
        switch nuclide
            case '10Be'
                Ps  = Ps10  + Ps10  .* truncnormrnd(1,Ps_uncert,-3*Ps_uncert,3*Ps_uncert);
                Pmu = Pmu10 + Pmu10 .* truncnormrnd(1,Pmu_uncert  ,-3*Pmu_uncert,3*Pmu_uncert);
                
                while any(Pmu10(:) < 0)   % if there are extreme values (negative ones) take a new sample
                    Pmu = Pmu10 + Pmu10 .* truncnormrnd(1,Pmu_uncert  ,-3*Pmu_uncert,3*Pmu_uncert);
                end
                
            case '36Cl'
                Ps   = Ps36  + Ps36    .* truncnormrnd(1,Ps_uncert,  -3*Ps_uncert  ,3*Ps_uncert);
                Pmu  = Pmu36 + Pmu36   .* truncnormrnd(1,Pmu_uncert, -3*Pmu_uncert ,3*Pmu_uncert);
                Pth  = Pth36  + Pth36  .* truncnormrnd(1,Pth_uncert, -3*Pth_uncert ,3*Pth_uncert);
                Peth = Peth36 + Peth36 .* truncnormrnd(1,Peth_uncert,-3*Peth_uncert,3*Peth_uncert); 
                
                while any(Pmu10(:) < 0 | Pth10(:) < 0 | Peth10(:) < 0)   % if there are extreme values (negative ones) take a new sample
                    Pmu = Pmu10 + Pmu10 .* truncnormrnd(1,Pmu_uncert  ,-3*Pmu_uncert,3*Pmu_uncert);
                    Pth  = Pth36  + Pth36  .* truncnormrnd(1,Pth_uncert, -3*Pth_uncert ,3*Pth_uncert);
                    Peth = Peth36 + Peth36 .* truncnormrnd(1,Peth_uncert,-3*Peth_uncert,3*Peth_uncert);
                end
                
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

            % Now that the number of layers for the section is known, an 
            % equal number of waiting times is drawn also from a power law 
            % distribution. To make sure that the total amount of time is
            % equal to the time covered by the section (i.e. the difference
            % between lower and upper age) the waiting times are rescaled to 
            % that duration
            % (!! need to check if this is allowed, i.e. is the distribution 
            % of the rescalled waiting times also a power law distribution?)
            dt = Intervals(j,2);
            pos_time = exp(exprnd(1/p_time,length(tmp_thickness),1));
            tmp_time = dt*pos_time/sum(pos_time);

            % The data of this section is attached on top of the data from the previous section
            Sed_col = [[tmp_thickness,tmp_time] ; Sed_col]; % bind all data together
        end
%         Sed_col = Sed_col(1:end-1,:);       % why remove the last????

        % Second part of the routine
        % The postdepositional nuclide production can be computed for the 
        % simulated core based on the depth/age constrain imposed by the 
        % Sed_col matrix.
        % The production of nuclides is computed from bottom to top (the 
        % loop indices decreases). Each part of the loop calculates the 
        % nuclide production of a given layer as well as all layers under it.
        % Each iteration is then summed with the previous one
        
        profile = zeros(1,section_depth);
        for k =  length(Sed_col):-1:1
            % For each layer temporary depth, density, time and concentration vectors are created
            % These vectors are made such that each element of the vector corresponds to 1 cm.
            tmp_depth = 1:sum(Sed_col(length(Sed_col):-1:k,1));
            tmp_time  = repmat(Sed_col(k,2), length(tmp_depth),1);
            
            % get mean depositional age of this layer for correction
            % scaling factor
            if k == length(Sed_col)
                tmp_age   = Depth_age_guess(end,2) - cumsum(Sed_col(k:length(Sed_col),2))/2;
            else
                tmp_age   = Depth_age_guess(end,2) - (sum(Sed_col(k:length(Sed_col),2))+sum(Sed_col(k+1:length(Sed_col),2)))/2;
            end
            
            Pind = round(tmp_age/1e2);  % get index of production rate column that is closest in age to tmp_age
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
            profile((section_depth-length(tmp_depth)+1):section_depth) = profile((section_depth-length(tmp_depth)+1):section_depth) + tmp_conc';
        end
        
        % add production after end of deposition
        if Depth_age_guess(1,2) ~= 0 % if there is a no deposition period at the end
            Pind = round(Depth_age_guess(1,2)/1e2);
            switch nuclide
                case '10Be'
                    profile(end) = profile(end) + (sum(Ps(section_depth,1:Pind) + ...
                        Pmu(section_depth,1:Pind)).*1e2);        % Production at/g
                case '36Cl'
                    profile(end) = profile(end) + (sum(Ps(section_depth,1:Pind) + ...
                        Pmu(section_depth,1:Pind) + Pth(section_depth,1:Pind) ...
                    + Peth(section_depth,1:Pind)).*1e2);        % Production at/g
            end
            profile(end) = profile(end).* exp(-1e2.*lambda);      % radioactive decay
        end
        
        PostProduction(i) = profile(end);
    % Increment progress bar
    waitbar(i/(n),wb)
end
close(wb)

Pmean = round(mean(PostProduction));
Pmedian = round(median(PostProduction));
Pstd = round(std(PostProduction));
Pquantiles = round(quantile(PostProduction,[0.17,0.83]));
disp(['postburial production ' name{1} ' = ' num2str(Pmean) ' +/- ' num2str(Pstd)]) 
disp(['postburial production median ' name{1} ' = ' num2str(Pmedian) ' +/- ' num2str(diff(Pquantiles)/2)]) 
disp(['postburial production quantiles ' name{1} ' = ' num2str(Pquantiles(1)) ' - ' num2str(Pquantiles(2))]) 


%% EXPORT THE DATA
if export
    filename = input('What name do you want to give this run for saving? ','s');
%     print(['PostburialProd_' name '_' filename '.eps'], '-depsc', '-painters') % save figure
%     savefig(h,['./output/PostburialProd_' ,name{1} ,'_' ,filename, '.fig'])

    % add all parameters to model
    Model.Ps_uncert   = Ps_uncert;
    Model.Pmu_uncert  = Pmu_uncert;
    Model.p_thickness = p_thickness;
    Model.p_time = p_time;
    Model.nruns = n;
    Model.mean_postburial = Pmean;
    Model.sd_postburial = Pstd;
    Model.median = Pmedian;
    Model.quants1783 = Pquantiles;

    save(['./output/PostburialProd_' name{1} '_' num2str(filename) '.mat'], 'Model')                % save model parameters
end
