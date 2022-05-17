% This is a function originally published in Marrero et al. 2016 and adapated for
% PostPro.
%%
%  age=cl36erateraw(pp,sp,sf,cp,maxerate,minrate)
%
%  Given the data for a saturated sample, computes the
%  corresponding erosion rate.
%
% This inner routine does not handle the uncertainty calculations,
% which are done by cl36erate.m.  Instead, this inner routine simply
% does the basic computation of the erosion rate.
%
% Inputs:
%    pp,sp,sf,cp       as from getpars36.
%    maxrate           Maximum erosion rate (g/(cm^2*kyr))
%    minrate           Minimum erosion rate (g/(cm^2*kyr))
%                      (optional, depfaults to 0)
%
% Returns
%
%   erate              g/(cm^2*kyr)
%
%
function erate=cl36erateraw(pp,sp,sf,cp,scaling_model,minrate)
%
% Set a default minimum erosion rate of 0 if none is specified.
%
if (nargin < 6)
  minrate=0.0;
end
%
% Figure out the maximum possible depth at which we'll ever need a
% production rate.
%
maxage=500;                       % 2000ka should be saturated for 36Cl.  By trial and error I found that 500 does only ~0.5% diff in the result but doubles the speed
%
sp.epsilon=minrate;
maxcon=predN36(pp,sp,sf,cp,maxage,scaling_model,1);
%
% Make sure that the concentration is reasonable.
%
if (sp.concentration36 > maxcon)
  warning('This sample is oversaturated based on min erosion rate');
  erate=NaN;
  return;
end
%
% The main loop does bisection search to find the corresponding
% erosion rate.
%
lowerrate=minrate;
%find a max erosion rate by starting at a low number, checking concentration, and
%increasing by an order of magnitude until the concentration is within the
%bounds of search.
MAXUPPERRATE=300; %this number is from trial and error; put in place to (80 in original)
% keep code from erroring out. Sept 2016 Shasta
upperrate=MAXUPPERRATE/1000; %divides it into 3 steps by order of magnitudesp.epsilon=upperrate;
minconc=predN36(pp,sp,sf,cp,maxage,scaling_model,1);

while (minconc-sp.concentration36 > 0)
    upperrate=upperrate*10;
    sp.epsilon=upperrate;
    %recalculate the maxdepth so that cp can be found to an appropriate depth
    maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+1000;
    %
    % Computed parameters.
    %
    cp=comppars36(pp,sp,sf,maxdepth);
    minconc=predN36(pp,sp,sf,cp,maxage,scaling_model,1);
    if upperrate > MAXUPPERRATE
       warning(['Maximum erosion rate (' num2str(MAXUPPERRATE) 'g/cm^2) reached - check sample inputs'])
       erate=NaN;
       return;
    end
end


while (upperrate-lowerrate > 1.0e-4)
  midrate=(upperrate+lowerrate)/2;
  sp.epsilon=midrate;
  midcon=predN36(pp,sp,sf,cp,maxage,scaling_model,1);
  if (midcon > sp.concentration36)
    lowerrate=midrate;
  else
    upperrate=midrate;
  end
end
erate=(upperrate+lowerrate)/2;
