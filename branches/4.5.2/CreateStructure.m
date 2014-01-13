%% Create a Structure object
function vars = CreateStructure(NumNodes, NumElements, NumMaterials,...
                                NumSections, Dimensions, NodesPerElement,... 
                                NumDOF,NumRestrainedDOF, NumSlavedDOF, NumGravityNodes,vars)

%% Set properties
vars.NumNodes = NumNodes;
vars.NumElements = NumElements;
vars.NumMaterials = NumMaterials;
vars.NumSections = NumSections;
vars.Dimensions = Dimensions;
vars.NodesPerElement = NodesPerElement;
vars.NumDOF = NumDOF;
vars.NumRestrainedDOF = NumRestrainedDOF;
vars.NumSlavedDOF = NumSlavedDOF;
vars.NumGravityNodes = NumGravityNodes;

% Calculate number of free degrees of freedom.
vars.NumFreeDOF = vars.NumDOF * vars.NumNodes-vars.NumRestrainedDOF-vars.NumSlavedDOF;
% vars.BOUND = zeros(NumNodes, NumDOF);
vars.EqDOF = zeros(NumNodes, NumDOF);

% initial condition for the time history analysis
vars.P0 = zeros(vars.NumFreeDOF,1);
