% This code calculates post burial production for 10Be
% and 36Cl (could easily expanded to work for 26Al and 14C). The
% calculations are based on geochronologic anchor ages. For sites with
% geochronologic ages the site specific parameters (latitude, elevation,
% shielding, if applicable chemistry etc) need to be specified. The code
% then takes the one or more samples and runs Monte Carlo simulation of
% deposition, where deposition occurs as linear aggradation and stops at a
% specifified age.
% Production rate profiles are calculated with CRONUCcalc v2.1 (Marrero et
% al. 2016)

% Richard Ott, 2021 (inspired by M. Lupker)
clc
clear
close all

addpath '.\subroutines'
addpath '..\..\..\..\Crete\Cretan_fans\data'

% USER CHOICE ----------------------------------------------------------- %
nuclide = '10Be';       % Choose '10Be' or '36Cl'
n = 1e3;                % number of runs
global scaling_model
scaling_model = 'lm';   % choose your scaling model, nomenclature follows Cronus

<<<<<<< HEAD
for i = 1:n
        % RANDOM SAMPLING ----------------------------------------------- %
        % Random sampling of the different parameters from a normal distribution for each of the n realisations of the simulation.
        % The ages of the core are sampled randomely within the uncertainty bounds of the dating provided to the function in Depth_age
%         Depth_age_guess = [normrnd(depth,depth_uncert),round(normrnd(age,age_uncert))];	% depth, age, matrix
        Depth_age_guess = [depth,round(normrnd(age,age_uncert))];	% depth, age, matrix
=======
% load sample data in Cronus excel format
% [num,txt,~] = xlsread('36Cl_data_CRONUS.xlsx','Matlab Postburial');
[num,txt,~] = xlsread('10Be_data_CRONUS','Matlab Postburial');
% load burial histories
[burial,burialtxt,~] = xlsread('Burialmodels','10Be_Model1');
% tag = 'test';
>>>>>>> linearAggradation_loop

inds = 1:2;   % indices of samples you intend to run (rows in input spreadsheet)
for i = inds
    %% assign data and constants ---------------1----------------------------- %

<<<<<<< HEAD
        % if there's an age inversion, take a new sample
        counter = 0;
        while any(diff(Depth_age_guess(:,2)) < 0)
%             Depth_age_guess = [normrnd(depth,depth_uncert),round(normrnd(age,age_uncert))];
            Depth_age_guess = [depth,round(normrnd(age,age_uncert))];
            counter = counter + 1;
            if counter > 1e4
                error("Couldn't sample a sequence without age inversion. Check yourage priors. ")
            end
        end
=======
    Perr = Puncerts(nuclide); % load uncertainties for production parameters
>>>>>>> linearAggradation_loop

    [Model,Para] = assignData(num,txt,burial,burialtxt,nuclide,inds(i)); % assign data
    if length(Para.depth) == 1 % skip modern samples that we not buried, remove for publication
        continue
    end

    %% PRODUCTION RATES ----------------------------------------------------- %

    Prod = ProductionParas(Para,nuclide);

    %% RUN FORWARD MODELS --------------------------------------------------- %

    Model = postburial_calc(Perr,Para,Model,Prod,nuclide,n);

    %% EXPORT THE DATA
%     save(['./output/PostburialProd_' Para.name{1} '_' tag '.mat'], 'Model')                % save model parameters
end
