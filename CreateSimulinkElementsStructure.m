function [ElementsStruct,ElementsStructBusName] = CreateSimulinkElementsStructure(Elements,Structure,RUNMODE)

ElementsStruct = struct; % create blank structure
% combine all elements into a giant structure
for element = 1:Structure.NumElements
    % Don't include Type 2 elements in the struct when doing an Experiment
    if Elements{element}.Type == 2 && strcmp(RUNMODE,'Experiment');
    else
        ElementsStruct.(sprintf('Element%d',element)) = Elements{element};
    end
end
% make a bus of this giant structure
ElementsStruct_bus = CreateBusWithDimensions(ElementsStruct);
ElementsStructBusName = ElementsStruct_bus.busName;


