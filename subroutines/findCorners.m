function [xMinP,xMaxP,yMinP,yMaxP] = findCorners(inGrid)
% the findCorners function takes an acsii grid input made with the function
% 'makeGrid' and finds the UTM coorinates of the corners of the grid.
% Created: 11/06/2013
% Author: Sean F. Gallen

ncols = inGrid.size(2);
nrows = inGrid.size(1);

[yMinP,xMinP] = pix2latlon(inGrid.refmat,nrows,1);
[yMaxP,xMaxP] = pix2latlon(inGrid.refmat,1,ncols);

end