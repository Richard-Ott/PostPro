function Perr = Puncerts(nuclide)

% Production rate uncertainties
switch nuclide
    case '10Be'
        Ps_uncert = 0.08;       % uncertainty spallation production 10Be (Phillips et al., 2016)
        Pmu_uncert = 0.3;       % uncertainty muon production 10Be (three times the muon uncertainty from Phillips 2016)
        Perr.Ps_uncert = Ps_uncert;
        Perr.Pmu_uncert = Pmu_uncert;
%     case '14C'
%         Ps_uncert = 0;          % uncertainty spallation production 14C (Phillips et al., 2016)
%         Pmu_uncert = 0.3;       % could not find a value for this yet
%         Perr.Ps_uncert = Ps_uncert;
%         Perr.Pmu_uncert = Pmu_uncert;
    case '36Cl'
        Ps_uncert  = 0.065;	    % relative uncertainty on the spallation production rate
        % average of PsCa and PsK (Cronus v2.1)
        Pmu_uncert = 0.25;      % relative uncertainty on the muon production rate
        % average between the f*Ca and f*K uncertainty in Marrero, 2016
        Pth_uncert  = 0.36; % relative uncertainty on the thermal neutron capture production rate (Cronus v2.1 with lm)
        Peth_uncert = 0.35; % relative uncertainty on the thermal neutron capture production rate (Cronus v2.1 with lm)
        
        Perr.Ps_uncert = Ps_uncert;
        Perr.Pmu_uncert = Pmu_uncert;
        Perr.Pth_uncert = Pth_uncert;
        Perr.Peth_uncert = Peth_uncert;
end
end

