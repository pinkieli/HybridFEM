%% Load Model Configuration
%  This reads the configuration and sets up the MATLAB workspace objects
%  Structure contains information about the overall structure
%  Elements contains the element information and its nodes and materials
%  Integrator contains properties of the integration method
[Structure, Elements, Integrator, Materials, Sections, Nodes, choice] = LoadConfiguration(INP_File,UICheck);
if strcmp(choice,'No')    
    fprintf('[FEM] Model Setup Aborted\n');
    break;
end

%% Calculate Matrices with respect to free degrees of freedom 
Structure = CalculateMatrices(Structure, Elements, Integrator);
NumFreeDOF = Structure.NumFreeDOF;  % Necessary for forming total RF

%% Eigen analysis to obtain structure periods
d=eig(Structure.StiffnessMatrixFree,Structure.MassMatrixFree);
Structure.System_periods = sort(2*pi*(ones(Structure.NumFreeDOF,1))./sqrt(d),'descend'); 
clear d;

%% Load the EQ Record
Integrator = GetEQHistory(EQ_File, Integrator, Structure, Elements);

%% Load Integration Parameters and set Initial Conditions
Integrator = InitializeIntegrator(Integrator, Structure);
Interpolations = Integrator.Interpolations;

%% Set up Timing and Data Storage
if strcmp(TARGET,'Simulink')
    %% Set sample rate and execution time for Simulink model
    SampleRate = Integrator.Timestep;   % Sample Rate of the Element blocks
    sample = SampleRate/Integrator.Interpolations;   % Sample rate of the Integrator/Controller
    RunningTime = SampleRate * (Integrator.Steps-1); % Set Running Time of the Simulink Model
    % Modify the EQ Record for Simulink compatibility
    Integrator.PEFFsimulink = horzcat([0:SampleRate:RunningTime-SampleRate]',Integrator.PEFF(2:end,:));     

    
else
    %% Initialize results vectors for Matlab simulation
    Am_Out_Steps = zeros(Integrator.Steps, 1);
    Am_Out_Displacement = zeros(Integrator.Steps, Structure.NumFreeDOF);
    Am_Out_RestoringForce = zeros(Integrator.Steps, Structure.NumFreeDOF);
    Am_Out_Velocity = zeros(Integrator.Steps,Structure.NumFreeDOF);
    Am_Out_Acceleration = zeros(Integrator.Steps,Structure.NumFreeDOF);    
    Am_Out_Acceleration(1,:) = Integrator.Acceleration;
end

%% Create and Initialize element recorder
ElementRecorderStruct = struct;
ElementRecorderStruct.Recorders = 0; % Need a dummy variable
if isnumeric(Structure.RecordElements) == 1
    Am_Out_ElementRestoringForce = CreateElementRecorder(Structure, Elements, Integrator);    
    for element = 1:length(Structure.RecordElements)
        elementType = Elements{Structure.RecordElements(element)}.Type;
        if (elementType == 9)
            ElementRecorderStruct.(sprintf('Element%d',Structure.RecordElements(element))) = zeros(12,1);
        else
            ElementRecorderStruct.(sprintf('Element%d',Structure.RecordElements(element))) = zeros(6,1);
        end
    end       
else    
    Am_Out_ElementRestoringForce = 0;    
end
ElementRecorderStruct_bus = CreateBusWithDimensions(ElementRecorderStruct);
ElementRecorderBus = eval(ElementRecorderStruct_bus.busName); 


%% Perform static analysis with gravity loadings
[Structure, Elements, Integrator, flag] = StaticAnalysis(Structure, Elements, Integrator);
if flag == 0    
    errordlg('Check Equilibrium under Static Loading','Analysis Error');
end
%% Store data to workspace
if strcmp(TARGET,'Matlab')
    % Store initial results vectors for Matlab Simulation
    Am_Out_Displacement(1,:) = Integrator.Displacement';
    Am_Out_RestoringForce(1,:) = Integrator.RestoringForce';
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
    [ElementsStruct,ElementsStructBusName] = CreateSimulinkElementsStructure(Elements,Structure,RUNMODE);
    ElementStructBus = eval(ElementsStructBusName);
    Elements = ElementsBack;
    clear ElementsBack;    
end


%% Print configuration to the screen
fprintf(['[FEM] Configuration   : ' INP_File,'\n']);
fprintf('[FEM] Structure\n');
fprintf('[FEM]  Nodes             : %i\n',Structure.NumNodes);
fprintf('[FEM]  Elements          : %i\n',Structure.NumElements);
fprintf('[FEM]  Materials         : %i\n',Structure.NumMaterials);
fprintf('[FEM]  Sections          : %i\n',Structure.NumSections);
fprintf('[FEM]  Restrained DOFs   : %i\n',Structure.NumRestrainedDOF);
fprintf('[FEM]  Slaved DOFs       : %i\n',Structure.NumSlavedDOF);
fprintf('[FEM]  Free DOFs         : %i\n',Structure.NumFreeDOF);
fprintf('[FEM]  Gravity Nodes     : %i\n',Structure.NumGravityNodes);
% Ask the User if they want to see structure periods?
if (UICheck == 1)
%     response = input('[FEM] View Structure Periods? y/[n]: ', 's');
    response = input('[FEM] View Mode shapes and Periods? y/[n]: ', 's');
else
    response = 'n';
end
if ~isempty(response) && ((response == 'y' || response == 'Y'))   
% The following commented out, since the mode shapes are being plotted (modified by CKolay)    
%     fprintf('[FEM] Structure Periods\n');
%     fprintf('[FEM]  %f sec\n',Structure.System_periods);
    NumModeShapes = input('Enter number of modes shapes to view = ');
    if NumModeShapes >0 && NumModeShapes <= Structure.NumFreeDOF
        PlotModeShapes(Structure, Elements, Nodes, NumModeShapes);
    else
        disp('Number of mode shapes cannot be negative or more than number of free DOF');
        break
    end
end
fprintf(['[FEM] EQ History      : ',EQ_File,'\n']);
fprintf('[FEM] EQ Scale Factor : %i\n',Integrator.EQScaleFactor);
fprintf('[FEM] Integrator      :');
switch (Integrator.MethodID)
    case 1
       fprintf(' CR\n'); 
    case 2
       fprintf(' Rosenbrock-W\n');
    case 3
       fprintf(' KR Semi-explicit alpha-method\n');
    case 4
       fprintf(' KR Explicit alpha-method\n');
    case 5
       fprintf(' Chang 2 Int Para\n');
    case 6
       fprintf(' Chang 3 Int Para\n');
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
    if (UICheck == 1)
        response = input('[FEM] Do you need to [u]pdate a Simulink Model file? u/[n]: ', 's');
    else
        response = 'n';
    end        
    % Update an existing model
    if ~isempty(response) && ((response == 'u' || response == 'U'))   
        % Create the element code
        if strcmp(RUNMODE, 'Simulation')
            [functionCode, elementCode] = SimulinkElementRFCodeGenerator(Elements, Structure, Integrator.MethodID, Am_Out_ElementRestoringForce, 1);
        end
        if strcmp(RUNMODE, 'Experiment')
            [functionCode, elementCode] = SimulinkElementRFCodeGenerator(Elements, Structure, Integrator.MethodID, Am_Out_ElementRestoringForce, 0);
        end 
        fprintf('[FEM] Updating existing model file\n');        
        [~, localfolder, ~] = fileparts(pwd);
        [filenameS, pathnameS] = uigetfile('*.mdl', 'Open existing HybridFEM Simulink Model File');
        if isequal(filenameS,0) || isequal(pathnameS,0)
            fprintf('[FEM] Model Update Canceled');
            return;
        else            
            existingModelFile = fullfile(pathnameS, filenameS);            
            ElementCodeInjector(existingModelFile, functionCode, elementCode);            
        end
    end
    % Ask the User if they want to export the xPC XML for the Element Recorder
    if (isnumeric(Structure.RecordElements) == 1 && length(Structure.RecordElements) > 0 && UICheck == 1)
        response = input('[FEM] Do you need to export the xPC XML for the Element Recorder? y/[n]: ', 's');
    else
        response = 'n';
    end     
    % Export the XML for the Element Recorder
    if ~isempty(response) && ((response == 'y' || response == 'Y'))  
        % Generate the code
        XMLforElementRecorder = '';
        for i = 1 : length(Structure.RecordElements)
            XMLforElementRecorder = [XMLforElementRecorder '<xPCSignal Location="0" Name="Element ' num2str(Structure.RecordElements(i)) ' N1UX" Units="n/a" Gain="1" isDAQ="False" />' char(10)];
            XMLforElementRecorder = [XMLforElementRecorder '<xPCSignal Location="0" Name="Element ' num2str(Structure.RecordElements(i)) ' N1UY" Units="n/a" Gain="1" isDAQ="False" />' char(10)];
            XMLforElementRecorder = [XMLforElementRecorder '<xPCSignal Location="0" Name="Element ' num2str(Structure.RecordElements(i)) ' N1UY" Units="n/a" Gain="1" isDAQ="False" />' char(10)];
            XMLforElementRecorder = [XMLforElementRecorder '<xPCSignal Location="0" Name="Element ' num2str(Structure.RecordElements(i)) ' N2UX" Units="n/a" Gain="1" isDAQ="False" />' char(10)];
            XMLforElementRecorder = [XMLforElementRecorder '<xPCSignal Location="0" Name="Element ' num2str(Structure.RecordElements(i)) ' N2UY" Units="n/a" Gain="1" isDAQ="False" />' char(10)];
            XMLforElementRecorder = [XMLforElementRecorder '<xPCSignal Location="0" Name="Element ' num2str(Structure.RecordElements(i)) ' N2UT" Units="n/a" Gain="1" isDAQ="False" />' char(10)];            
        end            
        [filenameS, pathnameS] = uiputfile('*.xml', 'Save to new XML file');
        if isequal(filenameS,0) || isequal(pathnameS,0)
            fprintf('[FEM] XML Generation Cancelled');
            return;
        else            
            xPCXMLforElementRecorderFile = fullfile(pathnameS, filenameS);                        
            fileID = fopen(xPCXMLforElementRecorderFile,'w');
            fprintf(fileID,'%s',char(XMLforElementRecorder));
            fclose(fileID);
            fprintf('[FEM] Element Recorder XML Generation Completed\n');
        end
    end
    fprintf('[FEM] Run Simulink model or compile for xPC now\n');
else
    input('[FEM] Hit Enter to begin numerical simulation');
    fprintf('[FEM] Running...\n');
    MatlabSimulation;
end