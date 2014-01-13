%% Create Element Type 4
function vars = CreateElement_Type4(ID,Type, Node1, Node2,...
    DampStiffFac, DampMassFac, Wi, m1, m2, mdof)

% Set element nodes
vars.ID = ID;
vars.Nodes = [Node1 Node2];
vars.Type=Type;
vars.DampStiffFac=DampStiffFac;
vars.DampMassFac=DampMassFac;

% Variables which are initialized from the input
vars.xyj(1,1) =  Node2.Xcoord;    
vars.xyj(2,1) =  Node2.Ycoord;
vars.xyi(1,1) =  Node1.Xcoord;
vars.xyi(2,1) =  Node1.Ycoord;

dx = vars.xyj(1,1) - vars.xyi(1,1); 
dy = vars.xyj(2,1) - vars.xyi(2,1);
L = sqrt(dx*dx+dy*dy);

vars.m1 = m1;
vars.m2 = m2;
vars.m = max([m1 m2])/1e10; % small mass 
vars.massDof = mdof;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change to dummy column for p-d effect
% for dummy column approach
% Wi is the weight assigned to node i and Hi is length of the column element 
if L < 1.0E-8
    vars.K22 = 0; vars.K25 = 0; vars.K52 = 0; vars.K55 = 0;
else
    vars.K22 =-Wi/L; vars.K25 = Wi/L; vars.K52 = Wi/L; vars.K55 =-Wi/L;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set handles to element matrices forming function and restoring force function
vars.FormElementMatrices = @FormElementMatrices_Type4;
vars.GetRestoringForce = @GetRestoringForce_Type4;