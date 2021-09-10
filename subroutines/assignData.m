function [Model,para] = assignData(num,txt,burial,burialtxt,nuclide)

ind  = input('Enter the number of the sample you want to run (the row within the CRONUS data file) ');

%% Load data                       
switch nuclide
    case '10Be'
        para.depth= num(ind,14)*100;        % depth of sample (cm), I enter the values as m into the excel sheet even though the Cronus input is in g/cmÂ²
        para.depth_uncert = num(ind,29)*100;% uncertainty in overburden (cm)
        para.elev = num(ind,3);             % elevation in (m)
        para.Tshd = num(ind,7);             % topographic shielding factor
        para.name = txt(ind+1,1);             % sample name
        para.age = num(ind,33)*1e3;         % age of sample   (yrs)
        para.age_uncert = num(ind,34)*1e3;  % age uncertainty (yrs)
        para.rho = num(ind,6);              % density of burying material n g/cm³

        % CONSTANTS ----------------------------------------------------- %
        para.lambda = log(2)/1.39e6;	       % decay constant for 10Be
        num(ind,4) = stdatm(para.elev);     % convert to air pressure
        
    case '36Cl'
        % CONSTANTS ----------------------------------------------------- %
        para.lambda = log(2)/3.013e5;       % decay constant for 36Cl (Audi, 2017)
        para.elev = num(ind,3);             % sample elvation
        para.name = txt(ind+1,1);             % sample name
        para.age = num(ind,80)*1e3;         % sample age in yrs
        para.age_uncert = num(ind,81)*1e3;
        para.depth = num(ind,12)*100;       % sample depth in cm
        para.depth_uncert = num(ind,51)*100;
        para.rho = num(ind,6);              % density of burying material in g/cm³

        % assign scaling factor for nucleonic and muonic production ----- %
        num(ind,4) = stdatm(para.elev);    % convert to air pressure
end

para.num = num(ind,:);  % just keep the sample you selected

% find the sample in the burial file
burialInds = find(strcmpi(para.name{1},burialtxt(:,1)))-1;  % find the indices of all layers within this stratigraphic section

% Add additonal constraints on deposition ------------------------------- %
% here we assemble the parameter of the layers that are not
% constrained by the sample you are running. This includes the surface and
% its age as well as additional layers in you section that may be dated and
% therefore constrain the deposition of overburden for your sample
para.depth         = [burial(burialInds,1)*100;para.depth];  % depths of constrained layers (cm), first should always be 0 for surface (cm)
para.depth_uncert = [burial(burialInds,2)*100;para.depth_uncert];  % (cm)
para.age           = [burial(burialInds,3)*1e3;para.age];  % ages of constrained layers in yrs
para.age_uncert   = [burial(burialInds,4)*1e3;para.age_uncert];  % age uncertainty of constrained layers in yrs

Model = struct();
Model.depths         = para.depth;          % depths of constrained layers (cm), first should always be 0 for surface (cm)
Model.depths_uncerts = para.depth_uncert;  % (cm)
Model.ages           = para.age;            % ages of constrained layers in yrs
Model.age_uncerts    = para.age_uncert;     % age uncertainty of constrained layers in yrs

end

