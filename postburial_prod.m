% This code calculates depth profiles of post burial production for 10Be
% and 36Cl (could easily expanded to work for 26Al and 14C). The
% calculations are based on geochronologic anchor ages. For the sites with
% geochronologic ages the site specific parameters (latitude, elevation,
% shielding, if applicable chemistry etc) need to be specified. The code
% then takes the one or more samples and runs Monte Carlo simulation of
% deposition, where deposition occurs as linear aggradation and stops at a
% specifified age.

% Richard Ott, 2021 (inspired by M. Lupker)
clc
clear
close all

addpath '.\subroutines'
addpath '..\..\..\..\Crete\Cretan_fans\data'

% USER CHOICE ----------------------------------------------------------- %
nuclide = '10Be';       % Choose '10Be' or '36Cl'
n = 5e3;                % number of runs
global scaling_model
scaling_model = 'lm';   % choose your scaling model, nomenclature follows Cronus

% load sample data in Cronus excel format
% [num,txt,~] = xlsread('36Cl_data_CRONUS.xlsx','Matlab Postburial');
[num,txt,~] = xlsread('10Be_data_CRONUS','Matlab Postburial');
% load burial histories
[burial,burialtxt,~] = xlsread('Burialmodels','10Be_model1');
tag = 'test';

inds = 1:7;
for i = inds
    %% assign data and constants ---------------1----------------------------- %

    Perr = Puncerts(nuclide); % load uncertainties for production parameters

    [Model,Para] = assignData(num,txt,burial,burialtxt,nuclide,inds(i)); % assign data

    %% PRODUCTION RATES ----------------------------------------------------- %

    Prod = ProductionParas(Para,nuclide);

    %% RUN FORWARD MODELS --------------------------------------------------- %

    Model = postburial_calc(Perr,Para,Model,Prod,nuclide,n);

    %% EXPORT THE DATA
    save(['./output/linearAggradation/' nuclide '/PostburialProd_' Para.name{1} '_' tag '.mat'], 'Model')                % save model parameters
end
