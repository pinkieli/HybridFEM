%% Get restoring force for the element
%  returns (6x1) restoring force vector of element in global coordinates
function [Element,RestoringForce] = GetRestoringForce_Type6(Element, Structure, Integrator)

% Get the displacement Ue (6x1) vector of the element in global coordinates
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
end

%ndof=6;  % The element has 6 dofs
du=Ue-Element.Uprev;  % The incremental displacement vector du of the element in global coordinates 
Element.Uprev=Ue; %Update Uprev of element for use in next state determination call

% Incremental basic deformations
dv = Element.avq * du;

% Integrate basic force and basic stiffness matrix
kb = zeros(3);
s  = zeros(3,1);

for i=1:length(Element.sections)
	
	% Form B matrix and compute section deformation
	B = [ Element.b1(i) Element.b2(i) 0; 0 0 Element.b3(i) ];
	dvs = B * dv;  
	
	% Section state determination and update section state
	[ss ks Element.sections(i)] = SectionState( Element.sections(i), dvs );
	
	% Sum section term
	kb = kb + B'*ks*B * ( Element.wi(i) * Element.L/2 );
	s  = s  + B'*ss   * ( Element.wi(i) * Element.L/2 );

end

% Restoring force and stiffness matrix
RestoringForce = Element.avq' * s;

if( Structure.RigidLinkNo )
    RestoringForce = TR' * RestoringForce;  
end
