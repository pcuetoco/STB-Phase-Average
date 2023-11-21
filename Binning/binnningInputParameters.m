function binningInputStruct = binnningInputParameters
%binnningInputParameters Return the input parameters needed to bin STB data
%   [binningInputStruct,phaseAverageInputStruct,pressureInputStruct] =
%   binnningInputParameters returns three structures containing the input
%   parameters for the binning of STB data. binningInputStruct contains the
%   standard time averaged binning input including binning parameters,
%   saving options project structure and visualization option.
%   phaseAverageInputStruct includes optional parameters for phase averaged
%   binning and pressureInputStruct includes optional parameters for
%   calculation of pressure from velocity data.

%% Standard time averaged binning input
%-------------------------- BINNING PARAMETERS --------------------------%
binningInputStruct.binSize = 18.76*1e-3;  % [m]
binningInputStruct.overlapFactor = 0.75;    % overlap factor
% Number of volume slices for binning parallelization
binningInputStruct.noSlices = 15;
% Ensemble averaging mode (see Aguera et al. 2016):
% - 'tophat'
% - 'gauss'
% - 'linear'
% - 'quadratic'
% - 'adaptive' --> for each bin it determines automatically the ensemble
% averaging mode based on the number of particles
binningInputStruct.averagingMethod = 'adaptive';
% Minimum number of particles for valid ensemble average
binningInputStruct.minNoParticles = 5;
% Option for cropping results to a volume of interest (comment if not
% desired) [xmin,xmax,ymin,ymax,zmin,zmax]
 binningInputStruct.domainCropVector = [-226.130,85.7416,-196.57,121.75,-150,150]*1e-3;  % [m]
%------------------------------------------------------------------------%
%-------------------------- SAVING OPTIONS --------------------------%
% Flag to store .mat files with x,u for each cone
%binningInputStruct.coneParticles2MatFilesFlag = false;
% Flag to store all particles' (x,u) in a single .mat file
%binningInputStruct.fullFieldParticles2MatFileFlag = false;
% Flag to store averaged velocity field in a single .mat file
%binningInputStruct.binnedData2MatFileFlag = false;
% Flag to write .dat file for tecplot (not recommended)
binningInputStruct.ascii4TecplotFlag = false;
% Flag to write .plt file for tecplot (recommended)
binningInputStruct.binary4TecplotFlag = true;
%--------------------------------------------------------------------%
%-------------------------- PROJECT STRUCTURE --------------------------%
% Path to folder of the project
binningInputStruct.projectFolderPath = ['C:\Users\pcuet\OneDrive\TU Delft\Thesis\05. W-Tunner First Tests\' ...
    '25042023_PropellerOnly\run_002'];
% Name of the folder containing the STB .dat files (one project may have
% different STB folders to distinguish between different calculations,
% leave the string empty if it STB folder coincides with project folder)
binningInputStruct.stbFolderName = '';
%-----------------------------------------------------------------------%
%------------------------- VISUALIZATION OPTIONS -------------------------%
binningInputStruct.plotParticlesFlag = false;
binningInputStruct.plotPlanesFlag = false;
%-------------------------------------------------------------------------%

%% Phase averaged binning input
% Temporal bin size defined as fraction of nominal period
binningInputStruct.temporalBinSize = 5*((1/(175.781/2))/360);   % [s]
% Number of desired phases
binningInputStruct.noPhases = 359;

%% Pressure input
% Pressure calculation flag: true = calculate pressure, false = don't
binningInputStruct.pressureFlag = false;
% Reference density [kg/m^3]
binningInputStruct.referenceDensity = 1.225;
% Reference freestream speed [m/s]
binningInputStruct.referenceFreestream = 11.8;
% Dirichlet flag: true = apply Dirichlet boundary condition, false = don't
binningInputStruct.dirichletFlag = true;
% Location(s) [x y z] of Dirichlet boundaries (array with different x,y,z
% along columns)
binningInputStruct.dirichletXyz = [-50,0,0]*1e-3;   % [m]
% Reference pressure for Dirichlet boundary condition
binningInputStruct.dirichletReferencePressure = 85;
end
