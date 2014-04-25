%% HybridFEM Simulation Script
%  This program loads a configuration file that contains information about
%  the structure, elements, nodes, materials and integration algorithm.
%  The analysis is completed through this script alone.

%% Set the files used for this program
%  INP_File: The Structure configuration file
%  EQ_File : The Earthquake ground acceleration record.  
%  Column 1 is time, Column 2 is the record.
%  This should be unscaled since the earthquake scaling factor is set in 
%  the INP_File. Otherwise, set the earthquake scaling factor to 1 in the 
%  INP_File. Gravity is not included in the effective force calculation.
%
%  TARGET can be either 'Matlab' or 'Simulink'.  The 'Matlab' option runs 
%  the analysis is completed through this script alone.  The 'Simulink'
%  option is used for preparing a RealtimeModel Simulink/xPC file.
%
%  RUNMODE can be either 'Simulation' or 'Experiment' and is only for 
%  when the TARGET is set to 'Simulink'.  This determines whether Elements
%  that are type 2 are using the inport for data from an  external source
%  such as a Load Cell or using a defined Stiffness value to calculate the 
%  Element's restoring force.
%
%  UICheck can be set to 0 or 1.  If it is set to 0, the user will not be
%  prompted to check parts of the setup.  
ClearWorkspace;
INP_File = '6StoryShearBldg_ExptSubsStructure_New.txt';
EQ_File = 'LOS270.txt';
TARGET = 'Simulink';
RUNMODE = 'Simulation';
UICheck = 1;  

%% Call Model Setup script
ModelSetup;

%%