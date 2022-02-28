This code package calculates postburial production for samples with complicated
time-depth burial histories. This code was developed for Ott et al. 2022 (see
citation below). Production rates are calculated using CRONUScalc v2.1 (Marrero
et al. 2016). The current version is developed for 10Be and 36Cl but can easily
be expanded to any nuclide within CRONUS.
To run the code you need to input your sample data through excel spreadsheets.
The nomenclature follows the CRONUScalc input for the nuclide samples. An
additional excel file with the parameters of the burial models (time-depth
constraints) needs to be provided.
To run the code CRONUScalc v2.1 needs to be in your Matlab path.

The postburial_prod.m script illustrates how to use the subroutines. The data
from Ott et al. 2022 are provided as example input data. 
Please report bugs to richard.ott1900@gmail.com

References:
Ott, Scherler etc 2022

Marrero, S.M., Phillips, F.M., Borchers, B., Lifton, N., Aumer, R., and Balco, G.
, 2016, Cosmogenic nuclide systematics and the CRONUScalc program: Quaternary
Geochronology, v. 31, p. 160–187, doi:10.1016/j.quageo.2015.09.005

LICENSE

Richard Ott, 2022
