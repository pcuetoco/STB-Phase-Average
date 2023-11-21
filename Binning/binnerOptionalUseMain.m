% Phase Average Binning
% 
% To run the code simply press either run buttons on the page after checking the parameters are set to the desired values.  The explanation of each of the parameters is at the end of the document.

% % Clear all variables, close all figure and files
clear
close all
fclose all;
clc

% Add Binning folder and subfolders to search path (substitute string with
% your local path)
currentFolder=pwd;
addpath(strcat(currentFolder,'\FileLoading'),strcat(currentFolder,'\Classes'),strcat(currentFolder,'\Utilities'));

%% Ask whether to load old .mat file for binning
answer = questdlg('Load old .mat file for binning? (Click No if .dat files needs to be read)',...
    'Load .mat file','Yes','No','No');
switch answer
    case 'Yes'
        % Select .mat file to load
        [matFilename,matFilepath] = uigetfile('*.mat',...
            'Select one .mat file');
        load(strcat(matFilepath,matFilename))
        % Retrieve binning input parameters
        %binningInputStruct = binnningInputParameters;

    case 'No'
 

        %% Select files to load
        % Select STB .dat files to be loaded

        stbDataPathArray = selectStbDatFiles(pwd);

           
        %% Load STB data and other files
        % Generate StbRun object 
        [particleXyzUvwArray,particleTimeArray] =...
            loadStbDataDavis10Light(stbDataPathArray);
        fprintf('Binned data imported. Proceeding with binning. ')

        %% Find reference initial time and cycle period for each STB run
        % Calculate offset of first timestep
        offset=particleTimeArray{1}(1,1); %[s]
        
        % Add time offset to current STB run
        particleTimeArray{1} = particleTimeArray{1}+(...
            -offset);
       
        

        %% Save workspace
        k = strfind(stbDataPathArray{1},filesep);
        saveFolderPath = stbDataPathArray{1}(1:k(end));
        binningInputStruct.saveFolderPath=saveFolderPath;
        filename = matlab.lang.makeValidName(['PreBinningLight',...
            datestr(now,'yyyymmdd')]);
        save([saveFolderPath,filename,'.mat'],'binningInputStruct',...
            'particleTimeArray','particleXyzUvwArray',...
            'saveFolderPath','-v7.3')
end
% 
% Enter binning parameters: 


        binningInputStruct.binSize =12*1e-3; %Bin size [m]
        binningInputStruct.overlapFactor = 0.75; %Overlap factor from 0 to 1

        % Number of volume slices for binning parallelization comment
        % optiona

        binningInputStruct.noSlices = 15;
        
        % averaging mode based on the number of particles
        binningInputStruct.averagingMethod = 'adaptive';
%  Optional cropping of the flowfield set all to zero to disable. 

          %  binningInputStruct.domainCropVector = [xmin,xmax,ymin,ymax,zmin,zmax]*1e-3;  % [m]
   

%% Phase averaged binning input
        % Temporal bin size defined as fraction of nominal period
        %Defined in degrees of freedom per phase.
        binningInputStruct.temporalBinSize = 5*((1/(175.781/2))/360);   % [s]
        % Number of desired phases
        binningInputStruct.noPhases = 30;

        % Minimum number of particles for valid ensemble average
        binningInputStruct.minNoParticles = 5;

        binningInputStruct.projectFolderPath= pwd;
        binningInputStruct.stbFolderName = '';
%Specify period
        period = 1/(175.781/2); %[s]
 % Calculate the nondimensional time vector of the current STB
        % run
        particleNondimensionalTimeArray{1,1} =...
            particleTimeArray{1}-floor(particleTimeArray{1}/...
            period)*period;

% Assemble final nondimensional time vector
        particleNondimensionalTimeVector = cell2mat(...
            particleNondimensionalTimeArray);

%% Perform binning over the different parameters defined
% Save the nondimensional time vector into binningInputStruct
binningInputStruct.particleNondimensionalTimeVector =...
    particleNondimensionalTimeVector;
% Set phase average flag to true
binningInputStruct.phaseAverageFlag = false;
binningInputStruct.useTemporalBinning = false;
binningInputStruct.zoneName="Test";

% 
%  Bin size: Size of the bin in mm. 
% OverlapFactor: Overlap factor between bins between 0 and 1
% Number of Slices:  Number of slices for pararell computing (Default is set to 15) beware, more slices does not neccesarilly mean it will compute faster as it might take longer to reconstruct. 
% Ensembe Averaging Mode: The ensemble averaging mode. Default is adaptive that changes the mode depending on available number of particles.  See : Agüera, N., Cafiero, G., Astarita, T., & Discetti, S. (2016). Ensemble 3D PTV for high resolution turbulent statistics. Measurement Science and Technology, 27(12), 124011. https://doi.org/10.1088/0957-0233/27/12/124011
% Optional Cropping: Set the desired size of the flowfield to compute. This feature is optional.
% Temporal Bin Size: Temporal bin size in seconds to use if temporal bins are enabled. 
% Number of Phases: Number of phases desired for output. 
% Min number of Particles: Minimum number of particles needed to consider a bin valid. 
% Project Folder Path: Location of the project where the results will be stored. The STB runs don't need to be in this same location.
% Phase Frequency: Frequency to perform the phase average in Hz.
% Phase Average Binning: Set to true for phase average, false for time average. 
% Temporal Binning Enabled: Set to true to use temporal bins, false if not desired. 
% Name to Display in Tecplot: The zone will have this name when oppened in tecplot. 





 
%% Perform Binning
binnedData = binningLight(particleXyzUvwArray,...
    binningInputStruct.binSize,binningInputStruct);
        
  
