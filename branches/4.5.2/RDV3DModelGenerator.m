%% HybridFEM RDV 3D Model Generation Script
%  A 3D RDV model is generated from the input file
%% Set the input file used for this program
%  INP_File: The Structure configuration file
if ~exist('INP_File','var')
    [INP_File, Pathname] = uigetfile('*.txt', 'Open the Input File');
    if isequal(INP_File,0) || isequal(Pathname,0)
        return;
    end
end

% Open 3D Model File
[MODEL_File, Pathname] = uiputfile('*.model', 'Save the Model File');
if isequal(INP_File,0) || isequal(Pathname,0)
    return;
else
    fid = fopen(MODEL_File,'w');
end
[~, MODEL_Name] = fileparts(MODEL_File);

%% Load Configuration
%  This reads the configuration and sets up the MATLAB workspace objects
%  Structure contains information about the overall structure
%  Elements contains the element information and its nodes and materials
[Model, Elements, Integrator, Materials, Sections, Nodes, choice] = LoadConfiguration(INP_File);

% Write Model File
fprintf(fid,'%s\n','# Begin Node section');
% Nodes
for i = 1:Structure.NumNodes
    fprintf(fid,'%i (%f, %f, %f) [%s,%s,%s]\n',i,Nodes(i).Xcoord,Nodes(i).Ycoord,Nodes(i).Zcoord,['Input/DOF' num2str(Nodes(i).UX)],['Input/DOF' num2str(Nodes(i).UY)],['Input/DOF' num2str(Nodes(i).THETA)]);
end
fprintf(fid,'%s\n','# End Node section');
fprintf(fid,'%s\n','===');
fprintf(fid,'%s\n','# Begin Member section');
% Elements
for i = 1:Structure.NumElements
    fprintf(fid,'%i %i\n',Elements{i}.Nodes(1).ID,Elements{i}.Nodes(2).ID);
end
fprintf(fid,'%s\n','# End Member section');
fprintf(fid,'%s\n','===');
fprintf(fid,'%s\n','# Begin ScaleNode section');
fprintf(fid,'%s\n','# End ScaleNode section');
fprintf(fid,'%s\n','===');
% End of Model Code

% Close File
fclose(fid);

fprintf('[FEM] 3D RDV Model Generation Completed\n');