function [OutRas1,OutRas2] = rescaleRasterModified(Raster1,Raster2,cropRaster,resRaster,method,shape)
% This function takes two rasters and makes them the same size and
% resolution. Input: cropRaster - number (1 or 2) to indicate which raster
% should be used as final size; rasRaster (1 or 2), indicate which raster
% should be the one for the final resolution. Method - the method you want
% for resolution resampling (e.g. 'bilinear','nearest'...).
% shape - true or false, do you want to nan all values thata re nan in the
% original DEM and therefore mimic its shape.
% This function requires the findcorners.m function!
% Richard Ott, 2019

% crop raster a bit wider and higher as supposed, to avoid boundary effects
% due to different cell sizes

safety = 1e4;           % added saftey distance for initial cropping
if cropRaster == 1
    [xMinP,xMaxP,yMinP,yMaxP] = findCorners(Raster1);
    Raster2 = crop(Raster2,[xMinP-safety,xMaxP+safety],[yMinP-safety,yMaxP+safety]);
elseif cropRaster == 2
    [xMinP,xMaxP,yMinP,yMaxP] = findCorners(Raster2);
    Raster1 = crop(Raster1,[xMinP-safety,xMaxP+safety],[yMinP-safety,yMaxP+safety]);
else
    disp('cropRaster must be numeric and either 1 or 2!')
end

% make them same resolution
if resRaster == 1
    Raster2 = resampleRas(Raster2,Raster1.cellsize,method);
elseif resRaster == 2
    Raster1 = resampleRas(Raster1,Raster2.cellsize,method);
else
    disp('resRaster must be numeric and either 1 or 2!')
end

% resize rasters
if cropRaster == 1
    Raster2 = crop(Raster2,[xMinP,xMaxP],[yMinP,yMaxP]);
elseif cropRaster == 2
    Raster1 = crop(Raster1,[xMinP,xMaxP],[yMinP,yMaxP]);
else
    disp('cropRaster must be numeric and either 1 or 2!')
end

% if desired make the outputs the same shape
if shape
    if cropRaster == 1
        Mask = isnan(Raster1.Z);
        Raster2.Z(Mask) = nan;
    else
        Mask = isnan(Raster2.Z);
        Raster1.Z(Mask) = nan;
    end
end

OutRas1 = Raster1;
OutRas2 = Raster2;
end

