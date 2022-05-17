% This is a function originally published in Marrero et al. 2016 and adapated for
% PostPro.
%
% function [Prodtotal,Prods,ProdsCa,ProdsK,ProdsTi,ProdsFe,Prodth,Prodeth,...
%    Prodmu,ProdmuCa,ProdmuK,Phith,Phieth,Kpercent,Capercent,Clpercent]=prodz36(z,pp,sf,cp);
%
% input:
%    z            A vector of depths (g/cm^2)
%    pp,sf,cp     Parameters for this sample, as returned by
%                 getpars36().
%
% The prodz function calculates the individual
% contribution to the production from spallation, epithermal,
% thermal, and muon pathways at a particular given depth.
%
%This code includes a snow shielding factor according to Zweck et al.
%(2013). This is set to 0 depth of cover normally, yielding a factor of 1.

function [Prodtotal,Prods,ProdsCa,ProdsK,ProdsTi,ProdsFe,Prodth,Prodeth,...
    Prodmu,ProdmuCa,ProdmuK,Phith,Phieth,Kpercent,Capercent,Clpercent]=prodz36(z,pp,sf,cp);
%
% Make sure that we don't have any depths that are too big.
%
if (max(z) > cp.maxdepth)
  disp('maxdepth is');
  cp.maxdepth
  disp('z is ');
  z
  error('Prodz called for depth greater than maxdepth.');
end
%
% Make sure that sf.currentsf.Sel36 (all of them) is a number.
%
if (isnan(sf.currentsf.Sel36Ca))
    error('sf.currentsf.Sel36Ca is a NaN');
end
if (isnan(sf.currentsf.Sel36K))
    error('sf.currentsf.Sel36K is a NaN');
end
if (isnan(sf.currentsf.Sel36Ti))
    error('sf.currentsf.Sel36Ti is a NaN');
end
if (isnan(sf.currentsf.Sel36Fe))
    error('sf.currentsf.Sel36Fe is a NaN');
end
if (isnan(sf.currentsf.SFth))
    error('sf.currentsf.SFth is a NaN');
end
if (isnan(sf.currentsf.SFeth))
    error('sf.currentsf.SFeth is a NaN');
end
%
% We'll get some slightly negative depths due to roundoff errors.
% Fix them.
%
for i=1:length(z)
  if (z(i)<-1.0e-4)
    error('z(i) is negative!');
  end
  if (z(i)<0.0)
    z(i)=0.0;
  end
end

%
%find the number of depths given in the input vector z
%
numberdepths=length(z);

%
%For each z, find the appropriate negative muon flux and total flux terms by
%interpolation
%assume depths are always given in an already-sorted vector
%

negfluxdepth=interpolate(cp.depthvector,cp.negflux,z);
totalfluxdepth=interpolate(cp.depthvector,cp.totalflux,z);
%
%create new variables for the depth-dependent muon/neutron variables used
%in thermal and epithermal production
%

Pmudepth=pp.Ys*(negfluxdepth)+0.0000058*(totalfluxdepth);
Rmudepth=Pmudepth/(sf.currentsf.SFeth*pp.Pf0*cp.Reth);
Rmupdepth=Rmudepth*pp.pEtha/cp.pEthss;

%
%depth independent variables for muon-produced neutrons
%

Pmu=pp.Ys*(cp.negflux(1))+0.0000058*(cp.totalflux(1));
Rmu=Pmu/(sf.currentsf.SFth*pp.Pf0*cp.Reth);
Rmup=Rmu*pp.pEtha/cp.pEthss;
%
%
% Snow shielding correction according to Zweck et al. (2013).
zsnow=0; %depth of snow cover [cm]
rhosnow=1; %density of snow [g/cm3]
zcover=zsnow*rhosnow; %snow mass length [g/cm2]
covertime=0.0; %fraction of time when cover is present. This is applied to
%each time step. Example: if a time step is 100 years, this would assume
%that the cover was present for half of that time.

%Zweck et al. (2013) constants for different composition rocks - leave only one option (a and
%b)uncommented. Uncertainties (sigma) - not used at this time. Assumes
%samples at surface.

% Siliceous dolomite
%a=[1.51 -0.428 0.37 740];
%b=[3.374 -0.0251 -2.228 -0.611 1.0166];
%sigmaa=[0.13 0.014 0.12 220];
%sigmab=[0.060 0.0013 0.069 0.040 0.0006];

% Basalt
%a=[1.87 -0.388 0.46 1000];
%b=[3.786 -0.0233 -2.604 -0.745 1.0194];
%sigmaa=[0.20 0.015 0.16 430];
%sigmab=[0.057 0.0011 0.073 0.046 0.0006];

% Granite
a=[1.81 -0.391 0.44 930];
b=[3.701 -0.0238 -2.525 -0.697 1.0185];
%sigmaa=[0.18 0.014 0.15 350];
%sigmab=[0.056 0.0011 0.069 0.041 0.0005];

snows=covertime*exp(-zcover/cp.Lambdafe)+(1-covertime); %Produces only 1 value - does
% not depend on depth or composition

snoweth=covertime*((a(1)*zcover+1)^(a(2))-(((cp.ls*cp.rb)*zcover^(a(3)))/a(4)))+(1-covertime);
% snoweth is depth-dependent so it produces a vector with one scaling
% factor for each depth

if zsnow>0
    snowth=covertime*((b(1)*exp(b(2)*zcover)+b(3)*exp(b(4)*zcover))*b(5).^(-(cp.ls*cp.rb)))+(1-covertime);
else
    snowth=1;
end
    % snowth is depth-dependent so it produces a vector with one scaling
% factor for each depth

%
% Introduce a new precomputed factor to speed up computation
%
expfactor=exp(-z/cp.Lambdafe);
%
% New productions for each spallation pathway separately
%
ProdsCa=sf.currentsf.Sel36Ca*sf.ST*snows*cp.PsCa*expfactor;
ProdsK=sf.currentsf.Sel36K*sf.ST*snows*cp.PsK*expfactor;
ProdsTi=sf.currentsf.Sel36Ti*sf.ST*snows*cp.PsTi*expfactor;
ProdsFe=sf.currentsf.Sel36Fe*sf.ST*snows*cp.PsFe*expfactor;

% Calculate total production via spallation
Prods=ProdsCa+ProdsK+ProdsTi+ProdsFe;
%
% Thermal neutron production
%
% First, get thermal neutron flux
%Old eqn: updated 1 Mar 2015 because the incoming flux is not the
%epithermal flux/thermal flux. Those are derived from the high-energy flux
%and so the scaling has to be the high-energy scaling. In this case, we use
%the general flux scaling factor over all energies.
%
Phith=sf.currentsf.SFth*sf.ST*snowth*(cp.Phistarthss*expfactor+...
     (1+Rmup).*cp.SFDeltaPhistarethss.*exp(-z/cp.Lethss)+...
     (1+Rmup.*cp.Rth).*cp.SFDeltaPhistarthss.*exp(-z/cp.Lthss)+...
     Rmupdepth.*cp.Phistarthss);
% Phith=sf.currentsf.SelSF*sf.ST*snowth*(cp.Phistarthss*expfactor+...
%      (1+Rmup).*cp.SFDeltaPhistarethss.*exp(-z/cp.Lethss)+...
%      (1+Rmup.*cp.Rth).*cp.SFDeltaPhistarthss.*exp(-z/cp.Lthss)+...
%      Rmupdepth.*cp.Phistarthss);
%
% Get thermal production using the flux
%
Prodth=cp.fth/cp.Lambdathss*Phith;
%
% Epithermal neutron production
%
% Calculate the flux of epithermal neutrons first
%Old eqn: updated 1 Mar 2015 because the incoming flux is not the
%epithermal flux/thermal flux. Those are derived from the high-energy flux
%and so the scaling has to be the high-energy scaling. In this case, we use
%the general flux scaling factor over all energies.

Phieth=sf.currentsf.SFeth*sf.ST*snoweth*(cp.Phistarethss.*expfactor+...
    (1+Rmu.*cp.Reth).*cp.FDeltaPhistareth.*exp(-z/cp.Lethss)+...
    Rmudepth.*cp.Phistarethss);
%Phieth=sf.currentsf.SelSF*sf.ST*snoweth*(cp.Phistarethss.*expfactor+...
%    (1+Rmu.*cp.Reth).*cp.FDeltaPhistareth.*exp(-z/cp.Lethss)+...
%    Rmudepth.*cp.Phistarethss);
%
% Get epithermal production using the flux
%
Prodeth=(1-cp.pEthss)*cp.feth/cp.Lambdaethss*Phieth;
%
% Muons are now being multiplied by the terrain shielding factor. They are
% already scaled to the site location.

ProdmuCa=sf.ST*interpolate(cp.muon36(1,:),cp.muon36(2,:),z);
ProdmuK=sf.ST*interpolate(cp.muon36(1,:),cp.muon36(3,:),z);
Prodmu=ProdmuCa+ProdmuK;

%
% Now, compute the total production.
%
Prodtotal=Prods+Prodth+Prodeth+Prodmu;
OneHundredOverProdtotal=100./Prodtotal;

%produce percent production for three main pathways
Kpercent=(ProdmuK+ProdsK).*OneHundredOverProdtotal;
Capercent=(ProdmuCa+ProdsCa).*OneHundredOverProdtotal;
Clpercent=(Prodth+Prodeth).*OneHundredOverProdtotal;



%A bit of debugging code to stop everything if NaN's have been
%produced anywhere in the computation.

if (sum(isnan(Prodtotal)) > 0)
  warning('Prodz36 produced NaN! Probably your erosion rate is high. production at greta depth will be set to zero. (RO edit)');
  Prodth(isnan(Prodth))   = 0;
  Prodeth(isnan(Prodeth)) = 0;
  Prodmu(isnan(Prodmu))   = 0;
  Prodtotal=Prods+Prodth+Prodeth+Prodmu;
end
