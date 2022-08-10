% 3D FDTD Code for MATLAB
% Author:   Stephen D. Gedney, University of Colorado Denver
%                              Electrical Engineering Department
% Date:     11/9/2015
% Distribution:   Authorized for distribution only to UCD students
% Disclaimer:     This code is not claimed to be error free.  Please
%                 report any bugs to stephen.gedney@ucdenver.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
tstart = tic;
% Pre-processing:
%--------------------------------------------
% Input mesh dimensions, materials, conformal boundaries, sources,
% time-signature, simulation control, output control, etc.
ReadInputData;

% Initializations:
%--------------------------------------------
% initialize the time-simulation (time-steps, simulation time)
InitilizeTimeSimulation;

% initialize main update coefficients:
InitializeUpdateCoefficients;

% pre-processing for subcell models such as conformal boundaries:
PreProcessSubCellGeometries;

% initialize PML parameters:
InitializePML;

% intialize the fields and auxiliary variables throughout the grid:
InitializeFields;

% intialize parameters for Near to far Field transformation
% Calculatro%%%ADDED
Initialize_FieldCalculator;


% TimeStepping:
%--------------------------------------------
% Advance the fields through nMax time iterations:
for n = 1:nMax
    % Explicit update of the electric fields:
    %---------------------------------------
    % update electric fields throughout the entire domain:
    MainEFieldUpdate;
    
    % Auxiliary electric field updates (material, local, PML):
    AuxiliaryEFieldUpdates;
    
    % Outer boundary update of E-fields (ABC, or if PEC do nothing)
    OuterBoundaryEFieldUpdates;
    
    % Electric field source injection:
    ETypeSourceInjection;
    
    % Explicit update of the magnetic fields:
    %---------------------------------------
    % update magnetic fields throughout the entire domain:
    MainHFieldUpdate;
    
    % Auxiliary magnetic field updates (material, local, PML):
    AuxiliaryHFieldUpdates;
    
    % Magnetic field source injection:
    HTypeSourceInjection;
    
    % Field calculators:
    %--------------------------------------
    % time-dependent field calculators (e.g., V, I, NF-FF)
    FieldCalculators;
    
    % Output Data
    OutputData;
    
    if(mod(n,100) == 0)
        fprintf('Completed %d out of %d time steps \n', n, nMax);
    end
end
% Post processing:
%--------------------------------------------
% perform all post-processing calculations
PostProcess;

tellapsed = toc(tstart);
fprintf('Simulation complete after %6.2f seconds \n', tellapsed');