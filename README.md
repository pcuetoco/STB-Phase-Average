# STB-Phase-Average
Matlab code to generate phase average and time-averaged binning based on STB.dat

## Usage

Download the all the files from the repository. The folder (~/ptv-binning-codes/Binning/) contains all the essential classes, functions and script for the ensemble averaging of your PTV data. In such folder, the script `phaseAverageBinnerMain.mlx` provides a script that you only need to fill in the parameters required for the binning and hit the run button. After running the program the interface will guide you through the steps for the binning. Beware to check the binning parameters before running the script. 

Alternatively, the script `binnerOptionalUseMain.m` can be used instead of the live script if preferred. It is recommended to use the livescript no to miss any required parameters, but the normal script is placed in case of need. 

The script will ask you to select the .dat file containing the particle information from the STB algorithm. Beware that reading a big STB file can take a while, and for this reason, the Matlab workspace is saved after reading the .dat file. After the first try it is possible to simply select the Matlab workspace, and the code will work with the new parameters that are set in the livescript. 

The resulting binning will be saved in a Tecplot binary file in the selected project folder. For each phase, a separate file will be saved. 

Don't hesitate to ask for questions if something is not correctly working. 


## Motivation
This repository contains the codes that I used to obtain ensemble averages of Particle Tracking Velocimetry (PTV) data for my [MSc Thesis] (https://repository.tudelft.nl/islandora/object/uuid:6d071306-faef-4cb9-9c21-a184fc8aaf0e?collection=education).

The original codes were developed for the following papers, please cite them if you use the codes of this repository to time-average your data:
- Jan F G Schneiders, Fulvio Scarano, Constantin Jux and Andrea Sciacchitano. "Coaxial volumetric velocimetry", 2018, [https://doi.org/10.1088/1361-6501/aab07d](https://doi.org/10.1088/1361-6501/aab07d)
- Constantin Jux, Andrea Sciacchitano, Jan F. G. Schneiders and Fulvio Scarano. "Robotic volumetric PIV of a full-scale cyclist", 2018, [https://doi.org/10.1088/1361-6501/aab07d](https://doi.org/10.1007/s00348-018-2524-1)

Or this one in case you phase-averaged your data: 

- Francesco M A Mitrotta, Jurij Sodja and Andrea Sciacchitano. "On the combined flow and structural measurements via robotic volumetric PTV", 2022, [https://doi.org/10.1088/1361-6501/ac41dd](https://doi.org/10.1088/1361-6501/ac41dd)
