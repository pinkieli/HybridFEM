%% Generate the code for the Element Calculation
function [functionCode, elementCode] = SimulinkElementRFCodeGenerator(Elements, Structure, MethodID, Out_ElementRestoringForce, isSimulation)

% Modify the function header for Element Recorder option
fileID = fopen('_functioncode.txt','w');
fprintf(fileID,'function [Total_RF ElementRFs] = eml(NumFreeDOF, Displacement, Velocity, ElementsStruct, RigidLinkNo, RigidLinkNodeID, RigidLinkMasterMatrix, RigidLinkSlaveMatrix)');    
% Close the temp file
fclose(fileID);
% Read it back to a string (faster)
functionCode = fileread([pwd '\_functioncode.txt']);
% Delete the temp file
delete([pwd '\_functioncode.txt']);
    

% Open a temp file to work with the element restoring force block
fileID = fopen('_elementcode.txt','w');

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
        % Record element restoring forces
        if isnumeric(Structure.RecordElements) == 1  
            for i = 1 : length(Structure.RecordElements)
                if Elements{element}.ID == Out_ElementRestoringForce(1,i).ID;
                    fprintf(fileID,'ElementRFs.Element%i = ElementRF;\\n',element);
                end
            end
        end       
    else
        if Elements{element}.Type == 2 && ~isSimulation ;
        else
            fprintf(fileID,'TotalRF = AssembleElementGlobalRF6DOF(Elements.Element%i, ElementRF, TotalRF);\\n',Elements{element}.ID);             
            % Record element restoring forces
            if isnumeric(Structure.RecordElements) == 1  
                for i = 1 : length(Structure.RecordElements)
                    if Elements{element}.ID == Out_ElementRestoringForce(1,i).ID;
                        fprintf(fileID,'ElementRFs.Element%i = ElementRF;\\n',element);
                    end
                end
            end       
        end
    end
end

% Close the temp file
fclose(fileID);
% Read it back to a string (faster)
elementCode = fileread([pwd '\_elementcode.txt']);
% Delete the temp file
delete([pwd '\_elementcode.txt']);