%% Create Element Type 8 
%Bouc Wen Model
function vars = CreateElement_Type8(ID, Type, Node1, Node2, DampStiffFac,...
                                    DampMassFac, dof, mat)
% Set element nodes
vars.ID = ID;
vars.Type=Type;           
vars.Nodes = [Node1 Node2];
vars.DampStiffFac=DampStiffFac;
vars.DampMassFac=DampMassFac;

% Set handles to element matrices forming function and restoring force
vars.FormElementMatrices = @FormElementMatrices_Type8;
vars.GetRestoringForce = @GetRestoringForce_Type8;
    
% keep track of element total deformation 
vars.Uprev=0.0;  

% Degrees of Freedom
vars.dof = dof;

% material model 
vars.prop = mat(1);      
