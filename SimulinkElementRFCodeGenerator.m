%% Generate the code for the Element Calculation
function elementCode = SimulinkElementRFCodeGenerator(Elements, Structure, isSimulation)

% Open a temp file
fileID = fopen('elementcode.txt','w');

for element = 1:Structure.NumElements
    % Only print out the Element type 2 if it is a Simulation
    if Elements{element}.Type == 2 && ~isSimulation ;    
    % Type 8 element uses Velocity
    elseif Elements{element}.Type == 8
            fprintf(fileID,'[Elements.Element%i, ElementRF] = GetRestoringForce_Type%i(Elements.Element%i, Displacement, Velocity, RigidLinkNo, RigidLinkNodeID, RigidLinkMasterMatrix, RigidLinkSlaveMatrix);\\n',Elements{element}.ID,Elements{element}.Type,Elements{element}.ID);
    % All other elements are normal
    else
        fprintf(fileID,'[Elements.Element%i, ElementRF] = GetRestoringForce_Type%i(Elements.Element%i, Displacement, RigidLinkNo, RigidLinkNodeID, RigidLinkMasterMatrix, RigidLinkSlaveMatrix);\\n',Elements{element}.ID,Elements{element}.Type,Elements{element}.ID);
    end
    
    % Print the Local to Global function
    if Elements{element}.Type == 9
        fprintf(fileID,'TotalRF = AssembleElementGlobalRF9DOF(Elements.Element%i, ElementRF, TotalRF);\\n',Elements{element}.ID);    
    else
        if Elements{element}.Type == 2 && ~isSimulation ;
        else
            fprintf(fileID,'TotalRF = AssembleElementGlobalRF6DOF(Elements.Element%i, ElementRF, TotalRF);\\n',Elements{element}.ID);
        end
    end
end

% Close the temp file
fclose(fileID);
% Read it back to a string (faster)
elementCode = fileread([pwd '\elementcode.txt']);
% Delete the temp file
delete([pwd '\elementcode.txt']);