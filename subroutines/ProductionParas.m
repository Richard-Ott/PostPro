function Prod = ProductionParas(para,nuclide)
% Use CRONUCScalc v2.1 Marrero et al. (2016) to calculate the production
% rate profiles with depth.
% Richard Ott, 2021

v2struct(para);
global scaling_model

pp=physpars();
max_depth = round(depth(end)*rho + 5*depth_uncert(end)*rho); % get maximum possible depth to come up in forward models, this will be the depth until which the production profile is calculated, in g/cm²
Prod.max_age   = age(end)   + 5*age_uncert(end);   % get maximum possible age for forward models, 5-sigma cutoff
switch nuclide
    case '10Be'
        num([14,29]) = 0;                          % set "depth-to-top-of-sample" for CRONUS back to zero in order to calculate full production rate profile
        [nominal10,uncerts10] = createage1026(num);               % get basic sample info 
        sp = samppars1026(nominal10);                             % Extract the sample parameters from the sampledatavector.
        sf = scalefacs1026(sp,scaling_model);                     % get scaling factors
        cp = comppars1026(pp,sp,sf,max_depth);                    % computed parameters
        
        % for all potential depths and ages calculate the production rate
        % profile. 
        % Production rates are sorted in matrix with depth rows (in g/cm2) and time
        % columns (every 100 years)
        Prod.Ps10  = nan(max_depth,ceil(Prod.max_age/1e2));
        Prod.Pmu10 = nan(max_depth,ceil(Prod.max_age/1e2));
        for i = 1:size(Prod.Ps10,2)
            % this scales back the productoin rate to a certain time, age has to be negative and in (ka)
            sf.currentsf=getcurrentsf(sf,-i*1e-1,scaling_model,'be');  
            % get production rates for a production profile
            [~,~,Prod.Ps10(:,i),Prod.Pmu10(:,i),~,~] = prodz1026(1:max_depth,pp,sf,cp);   
        end
        
    case '36Cl'
        num([12,51]) = 0;      % set "depth-to-top-of-sample" for CRONUS back to zero in order to calculate full production rate profile
        [nominal36,uncerts36,cov36]=createage36(num);    % Get basic info about the sample ages.
        sp = samppars36(nominal36);          % Extract the sample parameters from the sampledatavector.
        sf = scalefacs36(sp,scaling_model);  % get scaling factors
        cp = comppars36(pp,sp,sf,max_depth); % computed parameters
        
        % for all potential depths and ages calculate the production rate
        % profile. 
        % Production rates are sorted in matrix with depth rows (in g/cm2) and time
        % columns (every 100 years)
        Prod.Ps36   = nan(max_depth,ceil(Prod.max_age/1e2));
        Prod.Pmu36  = nan(max_depth,ceil(Prod.max_age/1e2));
        Prod.Pth36  = nan(max_depth,ceil(Prod.max_age/1e2));
        Prod.Peth36 = nan(max_depth,ceil(Prod.max_age/1e2));
        for i = 1:size(Prod.Ps36,2)
            % this scales back the production rate to a certain time, age has to be negative and in (ka)
            sf.currentsf=getcurrentsf(sf,-i*1e-1,scaling_model,'cl');     
            % get production rates for a production profile
            [~,Prod.Ps36(:,i),~,~,~,~,Prod.Pth36(:,i),Prod.Peth36(:,i),Prod.Pmu36(:,i),~,~,~,~,~,~,~] = ...
            prodz36(1:max_depth,pp,sf,cp);             
        end
end
end

