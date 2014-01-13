%% Get restoring force for the element
%  returns (6x1) restoring force vector of element in global coordinates
function [Element, RestoringForce] = GetRestoringForce_Type4(Element, Structure, Integrator)

% element geometry       
% element deformation % geometry restricted to 1,2 plane (x,y)
dx = Element.xyj(1,1) - Element.xyi(1,1); 
dy = Element.xyj(2,1) - Element.xyi(2,1);
L = sqrt(dx*dx+dy*dy);
dx = dx / L;
dy = dy / L;

% 6x6 transformation matrix from local to global coordiante
apq  = zeros(6);
apq(1,1) = dx; apq(2,2) = dx; apq(4,4) = dx; apq(5,5)= dx;
apq(1,2) = dy; apq(2,1) =-dy; apq(4,5) = dy; apq(5,4)=-dy;
apq(3,3) =1.0; apq(6,6) =1.0;

% element stiffness matrix in local coordinate
SE = zeros(6,6);  
SE(2, 2) = Element.K22;
SE(2, 5) = Element.K25;
SE(5, 2) = Element.K52;
SE(5, 5) = Element.K55;
% transformed element stiffness matrix in global coordinate
SE = apq' * SE * apq;

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
RestoringForce = SE*Ue;  % Resisting (6x1) force vector of element in global coordinates

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