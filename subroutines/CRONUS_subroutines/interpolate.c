/*
%
% yout=interpolate(xin,yin,x)
%
% Given an ordered set of (xin(i),yin(i)) pairs, with 
%   xin(1)< xin(2) < ... <xin(n)
%
% and an ordered set of x values
%
%   x(1) < x(2) < x(3) < ... < x(m)
%
% use linear interopolation to produce yout(1), yout(2), ..., yout(m).
%
% Returns NAN for any input x(i) where x(i) < xin(1) or
% x(i)>xin(n).
%
 */

#include <math.h>
#include "mex.h"
#define MAX(a,b) ((a)>(b)?(a):(b))

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  int i,j;
  int xinrows, xincolumns, xinlength, yinrows, yincolumns, yinlength;
  int xrows, xcolumns, xlength;  
  double *xin, *yin, *x, *yout;
  double NAN = strtod("NaN", NULL);

  if (nrhs != 3)
    mexErrMsgTxt("The number of input arguments must be 3.");

  if (nlhs != 1)
    mexErrMsgTxt("The number of output arguments must be 1.");

  xinrows=mxGetM(prhs[0]);
  xincolumns=mxGetN(prhs[0]);

  yinrows=mxGetM(prhs[1]);
  yincolumns=mxGetN(prhs[1]);

  xinlength=MAX(xinrows,xincolumns);
  yinlength=MAX(yinrows,yincolumns);

  if (xinlength != yinlength)
    mexErrMsgTxt("xin and yin must be of same length.");

  xrows=mxGetM(prhs[2]);
  xcolumns=mxGetN(prhs[2]);

  xlength=MAX(xrows,xcolumns);

  plhs[0] = mxCreateDoubleMatrix(xrows, xcolumns, mxREAL);

  /*
   * Get pointers to the actual data.
   */

  xin=mxGetPr(prhs[0]);
  yin=mxGetPr(prhs[1]);
  x=mxGetPr(prhs[2]);
  yout=mxGetPr(plhs[0]);

  /*
   * Now do the actual interpolate.
   */

  i=0;
  j=0;
  while (i<xlength)
    {
      if (x[i] < xin[0])
	{
	  yout[i]=NAN;
	  i++;
	}
      else
	{
	  if (x[i]>xin[xinlength-1])
	    {
	      yout[i]=NAN;
	      i++;
	    }
	  else
	    {
	      while (xin[j+1]<x[i])
		j++;
	      yout[i]=yin[j]+(yin[j+1]-yin[j])*(x[i]-xin[j])/(xin[j+1]-xin[j]);
	      i++;
	    };
	};
    };
}
