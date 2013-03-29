%% Create Element Type 2
function vars = CreateElement_Type2(ID,Type, Node1, Node2, DampStiffFac, DampMassFac, Stiffness, Cexp)

% Set element nodes
vars.ID = ID;
vars.Type=Type;
vars.Nodes = [Node1 Node2];
vars.DampStiffFac=DampStiffFac;
vars.DampMassFac=DampMassFac;

% Set element properties
vars.Stiffness = Stiffness;
vars.Cexp=Cexp;

% Set handles to element matrices forming function and restoring force
vars.FormElementMatrices = @FormElementMatrices_Type2;
vars.GetRestoringForce = @GetRestoringForce_Type2;