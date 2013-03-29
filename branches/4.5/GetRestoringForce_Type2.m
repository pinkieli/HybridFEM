%% Get restoring force for the element
%  returns (6x1) restoring force vector of element in global coordinates
function [Element, RestoringForce] = GetRestoringForce_Type2(Element, Structure, Integrator)

% Preset arrays
SE = zeros(6,6);
SE(1, 1) = Element.Stiffness;
SE(1, 4) = -Element.Stiffness;
SE(4, 4) = Element.Stiffness;
SE(4, 1) = SE(1,4);

% Create the displacement Ue (6x1) vector of the element in global coordinates
Ue = zeros(6,1);
dof=zeros(6,1);
dof(1,1)=Element.Nodes(1).UX;
dof(2,1)=Element.Nodes(1).UY;
dof(3,1)=Element.Nodes(1).THETA;
dof(4,1)=Element.Nodes(2).UX;
dof(5,1)=Element.Nodes(2).UY;
dof(6,1)=Element.Nodes(2).THETA;
for i=1:6
    if dof(i,1) ~= -1
        Ue(i,1)=Integrator.Displacement(dof(i,1),1);
    end
end
RestoringForce = SE*Ue;  % Resisting (6x1) resisting force vector of element in global coordinates

if( Structure.RigidLinkNo )
    TR = eye(6,6); % Transformation matrix for rigid link
    for jjj=1:2
        idx=find(Structure.RigidLinkNodeID(:,2) == Element.Nodes(jjj).ID); % find if a slave node exists in the element
        if ( idx )
            TL=CalculateRigidLinkTransformation(Structure.RigidLinkMaster{idx},Structure.RigidLinkSlave{idx} );
            TR(3*(jjj-1)+1:3*jjj,3*(jjj-1)+1:3*jjj) = TL;
        end
    end 
    
    Ue = TR * Ue;     
    RestoringForce = TR' * SE * Ue;  
end