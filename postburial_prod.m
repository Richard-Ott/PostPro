clc
clear
close all

addpath '.\subroutines'

% USER CHOICE ----------------------------------------------------------- %
nuclide = '10Be';       % Choose '10Be' or '36Cl'
n =5e2;                % number of runs
global scaling_model
scaling_model = 'lm';   % choose your scaling model, nomenclature follows Cronus

% [num,txt,~] = xlsread('36Cl_data_CRONUS.xlsx','Matlab Postburial');% load sample data in Cronus excel format
[num,txt,~] = xlsread('10Be_data_CRONUS','Matlab Postburial');       % load sample data in Cronus excel format
[burial,burialtxt,~] = xlsread('Burialmodels','10Be_Model1');        % load burial histories

%% LOAD DATA
Perr = Puncerts(nuclide); % load uncertainties for production parameters

[Model,Para] = assignData(num,txt,burial,burialtxt,nuclide,5); % assign data


%% PRODUCTION RATES ----------------------------------------------------- %

Prod = ProductionParas(Para,nuclide);

%% RUN FORWARD MODELS --------------------------------------------------- %

Model = postburial_calc(Perr,Para,Model,Prod,nuclide,n,'plot');


