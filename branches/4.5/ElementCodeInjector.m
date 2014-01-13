%% Replace the Simulink Element Code in a model file
function ElementCodeInjector(modelFile, functionCode, elementCode)
%% Element Function Code
% Declare the start and end strings to look for in the Element Block to 
%  replace the element functions
startString = java.lang.String('%*fs');
endString = java.lang.String('%*fe');

% read in mdl file text
text = fileread(modelFile);
modelFileData = java.lang.String(text);

% find index of start, save beginning
startIndex = modelFileData.indexOf(startString);
if startIndex < 1
    fprintf('[FEM] Error: Not a HybridFEM Simulink File\n');
    return;
end
firstHalf = modelFileData.substring(0, startIndex+startString.length()+2);

% find index of end, save end
endIndex = modelFileData.indexOf(endString);
if endIndex < 1
    fprintf('[FEM] Error: Not a HybridFEM Simulink File\n');
    return;
end
lastHalf = modelFileData.substring(endIndex-1);

% combine beginning, new code and end
if isempty(functionCode)   
    fprintf('[FEM] Error: No Element Function Code Supplied\n');
    return;
end
newModelFile = firstHalf.concat(functionCode);
newModelFile = newModelFile.concat(lastHalf);

%% Element Block Code
% Declare the start and end strings to look for in the Element Block to 
%  replace the element functions
startString = java.lang.String('%*s');
endString = java.lang.String('%*e');

% read in mdl file text
%text = fileread(modelFile);
%modelFileData = java.lang.String(text);

% find index of start, save beginning
%startIndex = modelFileData.indexOf(startString);
startIndex = newModelFile.indexOf(startString);
if startIndex < 1
    fprintf('[FEM] Error: Not a HybridFEM Simulink File\n');
    return;
end
%firstHalf = modelFileData.substring(0, startIndex+startString.length()+2);
firstHalf = newModelFile.substring(0, startIndex+startString.length()+2);

% find index of end, save end
%endIndex = modelFileData.indexOf(endString);
endIndex = newModelFile.indexOf(endString);
if endIndex < 1
    fprintf('[FEM] Error: Not a HybridFEM Simulink File\n');
    return;
end
%lastHalf = modelFileData.substring(endIndex-2);
lastHalf = newModelFile.substring(endIndex-2);

% combine beginning, new code and end
if isempty(elementCode)   
    fprintf('[FEM] Error: No Element Code Supplied\n');
    return;
end
newModelFile2 = firstHalf.concat(elementCode);
newModelFile2 = newModelFile2.concat(lastHalf);

%% Save new model to file
[filename, pathname] = uiputfile('*.mdl', 'Save new MDL file');
if isequal(filename,0) || isequal(pathname,0)
    fprintf('[FEM] New Model Generation Canceled\n');
    return;
else
    newFile = fullfile(pathname, filename);
    fileID = fopen(newFile,'w');
    fprintf(fileID,'%s',char(newModelFile2));
    fclose(fileID);
    fprintf('[FEM] New Model Generation Completed\n');
end
