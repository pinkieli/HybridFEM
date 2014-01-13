%% Load Model Configuration
%  This reads the configuration and sets up the MATLAB workspace objects
%  Structure contains information about the overall structure
%  Elements contains the element information and its nodes and materials
%  Integrator contains properties of the integration method
[Structure, Elements, Integrator, Materials, Sections, Nodes, choice] = LoadConfiguration(INP_File);
if strcmp(choice,'No')    
     fprintf('[FEM] Model Setup Aborted\n');
     break;
end

%% Calculate Matrices with respect to free degrees of freedom 
Structure = CalculateMatrices(Structure, Elements, Integrator);
NumFreeDOF = Structure.NumFreeDOF;  % Necessary for forming total RF

%% Eigen analysis to obtain structure periods
d=eig(Structure.StiffnessMatrixFree,Structure.MassMatrixFree);
Structure.System_periods = 2*pi*(ones(Structure.NumFreeDOF,1))./sqrt(d); 
clear d;

%% Load the EQ Record
Integrator = GetEQHistory(EQ_File, Integrator, Structure, Elements);

%% Load Integration Parameters and set Initial Conditions
Integrator = InitializeIntegrator(Integrator, Structure);
Interpolations = Integrator.Interpolations;

%% Set up Timing and Data Storage
if strcmp(TARGET,'Simulink')
    %% Set sample rates and execution time for Simulink model
    SampleRate = Integrator.Timestep;   % Sample Rate of the Element blocks
    sample = SampleRate/Integrator.Interpolations;   % Sample rate of the Integrator/Controller
    RunningTime = SampleRate * Integrator.Steps; % Set Running Time of the Simulink Model
    % Modify the EQ Record for Simulink compatibility
    %Integrator.PEFFrt = horzcat([SampleRate:SampleRate:RunningTime]',Integrator.PEFF([2:end],:));     
    Integrator.PEFFsimulink = horzcat([SampleRate:SampleRate:RunningTime]',Integrator.PEFF);     
else
    %% Initialize results vectors for Matlab simulation
    Out_Steps = zeros(Integrator.Steps, 1);
    Out_Displacement = zeros(Integrator.Steps, Structure.NumFreeDOF);
    Out_RestoringForce = zeros(Integrator.Steps, Structure.NumFreeDOF);
    Out_Velocity = zeros(Integrator.Steps,Structure.NumFreeDOF);
    Out_Acceleration = zeros(Integrator.Steps,Structure.NumFreeDOF);    
end

%% Perform static analysis with gravity loadings
[Structure, Elements, Integrator, flag] = StaticAnalysis(Structure, Elements, Integrator);
if flag == 0    
    errordlg('Check Equilibrium Under static loading','Analysis Error');
end

%% Store data to workspace
if strcmp(TARGET,'Matlab')
    %% Store initial results vectors for Matlab Simulation
    %Out_Displacement(1,:) = Integrator.Displacement';
    %Out_RestoringForce(1,:) = Integrator.RestoringForce';
else     
    % Save RigidLink data to workspace for Simulink
    RigidLinkNo = Structure.RigidLinkNo;
    RigidLinkNodeID = Structure.RigidLinkNodeID;
    RigidLinkMasterMatrix = Structure.RigidLinkMasterMatrix;
    RigidLinkSlaveMatrix = Structure.RigidLinkSlaveMatrix;
    
    % Determine Run Mode for Simulink Element 2
    if strcmp(RUNMODE,'Experiment')
        RUNMODE_BIT = 1;
    else
        RUNMODE_BIT = 0;
    end
end

%% Create Element and Section Arrays
if strcmp(TARGET,'Simulink')
    % Save a copy of the Elements before the generic Simulink Bus is generated
    ElementsBack = Elements;
    % Remove references to functions for compatibility with Simulink
    for element = 1:Structure.NumElements
        Elements{element} = rmfield(Elements{element},'FormElementMatrices');
        Elements{element} = rmfield(Elements{element},'GetRestoringForce');        
    end
%    ElementsArray = CreateSimulinkElementArray(Elements, Structure);    
    [ElementsStruct,ElementsStructBusName] = CreateSimulinkElementsStructure(Elements,Structure,RUNMODE);
    ElementStructBus = eval(ElementsStructBusName);
    Elements = ElementsBack;
    clear ElementsBack;
end

%% Print configuration to the screen
fprintf(['[FEM] Configuration   : ' INP_File,'\n']);
fprintf('[FEM] Structure\n');
fprintf('[FEM]  Nodes          : %i\n',Structure.NumNodes);
fprintf('[FEM]  Elements       : %i\n',Structure.NumElements);
fprintf('[FEM]  NumFreeDOF     : %i\n',Structure.NumFreeDOF);
% Ask the User if they want to see structure periods?
response = input('[FEM] View Structure Periods? y/[n]: ', 's');
if ~isempty(response) && ((response == 'y' || response == 'Y'))   
    fprintf('[FEM] Structure Periods\n');
    fprintf('[FEM]  %f sec\n',Structure.System_periods); 
end
fprintf(['[FEM] EQ History      : ',EQ_File,'\n']);
fprintf('[FEM] EQ Scale Factor : %i\n',Integrator.EQScaleFactor);
fprintf('[FEM] Integrator      :');
switch (Integrator.MethodID)
    case 1
       fprintf(' CR\n'); 
    case 2
       fprintf(' Rosenbrock-W\n');
end
fprintf('[FEM] Steps           : %i\n',Integrator.Steps);
fprintf('[FEM] Timestep        : %f seconds\n',Integrator.Timestep);
fprintf('[FEM] Target          : %s \n',TARGET);
if strcmp(TARGET,'Simulink')
    if strcmp(RUNMODE, 'Simulation')
        fprintf('[FEM] Run Mode        : Simulation\n');
    end
    if strcmp(RUNMODE, 'Experiment')
        fprintf('[FEM] Run Mode        : Experiment\n');
    end    
    fprintf('[FEM] Sample Rate     : %f seconds\n',sample);
    fprintf('[FEM] Running Time    : %f seconds\n',RunningTime);    
    % Ask the User if they want to Create a new model
    response = input('[FEM] Do you need to [c]reate or [u]pdate a Simulink Model file? c/u/[n]: ', 's');
    % Create a new model   
    % Create the element code
    if strcmp(RUNMODE, 'Simulation')
        elementCode = SimulinkElementRFCodeGenerator(Elements, Structure, 1);
    end
    if strcmp(RUNMODE, 'Experiment')
        elementCode = SimulinkElementRFCodeGenerator(Elements, Structure, 0);
    end 
    if ~isempty(response) && ((response == 'c' || response == 'C'))   
        fprintf('[FEM] Copying model file from template\n');        
        [~, localfolder, ~] = fileparts(pwd);
        [filenameS, pathnameS] = uigetfile('*.mdl', 'Locate SimulinkModel_Template.mdl in HybridFEM folder');
        if isequal(filenameS,0) || isequal(pathnameS,0)
            fprintf('[FEM] New Model Generation Canceled\n');
            return;
        else            
            templateModelFile = fullfile(pathnameS, filenameS);            
            ElementCodeInjector(templateModelFile, elementCode);            
        end
    end
    % Update an existing model
    if ~isempty(response) && ((response == 'u' || response == 'U'))   
        fprintf('[FEM] Updating existing model file\n');        
        [~, localfolder, ~] = fileparts(pwd);
        [filenameS, pathnameS] = uigetfile('*.mdl', 'Open existing HybridFEM Simulink Model File');
        if isequal(filenameS,0) || isequal(pathnameS,0)
            fprintf('[FEM] Model Update Canceled');
            return;
        else            
            existingModelFile = fullfile(pathnameS, filenameS);            
            ElementCodeInjector(existingModelFile, elementCode);            
        end
    end
    fprintf('[FEM] Run Simulink model or compile for xPC now\n');
else
    input('[FEM] Hit Enter to begin numerical simulation');
    fprintf('[FEM] Running...\n');
    MatlabSimulation;
end