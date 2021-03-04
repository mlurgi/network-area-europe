# network-area-europe

Start here!

Authors: Dr Miguel Lurgi
Lecturer in Biosciences (Computational Ecology)
Computational Ecology Lab - Department of Biosciences
Swansea University, UK
 
and

Centre for Biodiversity Theory and Modelling
Theoretical and Experimental Ecology Station, CNRS, France

AND

Dr Núria Galiana
Centre for Biodiversity Theory and Modelling
Theoretical and Experimental Ecology Station, CNRS, France

Date Created: 19-12-2020

Copyright (c) Miguel Lurgi, 2020
Email: miguel.lurgi@swansea.ac.uk

This repository contains the source code used to produce the results presented in our paper entitled:

'The spatial scaling of food web structure across European biogeographical regions'

Published in Ecography on February 2021 - https://doi.org/10.1111/ecog.05229

In this file I describe the contents of the repository and provide a brief description of what source file does when executed.

To execute the scripts successfully, place all the scripts, folders, and data files contained in this repository within the same folder.

The main script file is:

1.- 'european-network-area-relationships.r' : This script implements different map cells aggregation algorithms to simulate increases in spatial scales across Europe and calculates food web properties on the food web determined by the species present in each given area. To do this it manipulates species distribution maps and a meta-network (metaweb) of trophic interactions constructed at the European regional scale.

In addition to the main script, another source code file containing helper functions is provided to support calculations of food web properties: utils.r. A data file is also provided to facilitate running of the script, test-data.rda, which contains a sample of the data used for the analysis. The folders BiogeoRegions2016_shapefile and mask10k are also provided for completeness. These are needed to retrieve the spatial data necessary to build communities and networks across the European geographical extent.

I hope you enjoy the code!

If you have any questions / run into any issues, please contact me at miguel.lurgi@swansea.ac.uk





