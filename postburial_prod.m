clc
clear
close all

addpath '.\subroutines'
addpath '..\..\..\..\Crete\Cretan_fans\data'

% USER CHOICE ----------------------------------------------------------- %
nuclide = '36Cl';       % Choose '10Be' or '36Cl'
n = 1e2;                % number of runs
global scaling_model
scaling_model = 'lm';   % choose your scaling model, nomenclature follows Cronus

[num,txt,~] = xlsread('36Cl_data_CRONUS.xlsx','Matlab Postburial');% load sample data in Cronus excel format
% [num,txt,~] = xlsread('10Be_data_CRONUS','Matlab Postburial');       % load sample data in Cronus excel format
[burial,burialtxt,~] = xlsread('Burialmodels','36Cl_Model1');        % load burial histories
% tag = 'test';

inds = 4;   % indices of samples you intend to run (rows in input spreadsheet)
for i = 1:length(inds)
    %% assign data and constants ---------------1----------------------------- %

    Perr = Puncerts(nuclide); % load uncertainties for production parameters

    [Model,Para] = assignData(num,txt,burial,burialtxt,nuclide,inds(i)); % assign data
    if length(Para.depth) == 1 % skip modern samples that we not buried, can be removed in most cases
        continue
    end

    %% PRODUCTION RATES ----------------------------------------------------- %

    Prod = ProductionParas(Para,nuclide);

    %% RUN FORWARD MODELS --------------------------------------------------- %

    Model = postburial_calc(Perr,Para,Model,Prod,nuclide,n,'plot');

    %% EXPORT THE DATA
%     save([Para.name{1} '_' tag '.mat'], 'Model')                % save model parameters
end
