%% Get restoring force for the element
%  returns (6x1) restoring force vector of element in global coordinates
function [Element,RestoringForce] = GetRestoringForce_Type8(Element,Structure,Integrator)
% Get the displacement Ue (6x1) vector of the element in global coordinates
Uev = zeros(6,1);
Vev = zeros(6,1);
dof=zeros(6,1);
dof(1,1)=Element.Nodes(1).UX;
dof(2,1)=Element.Nodes(1).UY;
dof(3,1)=Element.Nodes(1).THETA;
dof(4,1)=Element.Nodes(2).UX;
dof(5,1)=Element.Nodes(2).UY;
dof(6,1)=Element.Nodes(2).THETA;
for i=1:6
    if dof(i,1) ~= -1
        Uev(i,1)=Integrator.Displacement(dof(i,1),1);
        Vev(i,1)=Integrator.Velocity(dof(i,1),1);        
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
    Uev = TR * Uev;         
    Vev = TR * Vev;
end

% deformation and velocity
% if Element.dof ==1
    Ue = Uev(Element.dof,1)-Uev(Element.dof+3,1);
    Ve = Vev(Element.dof,1)-Vev(Element.dof+3,1);
% elseif Element.dof ==2
%     Ue = Uev(2,1)-Uev(5,1);
%     Ve = Vev(2,1)-Vev(5,1);
% else
%     Ue = Uev(3,1)-Uev(6,1);
%     Ve = Vev(3,1)-Vev(6,1);
% end

% Incremental displacement and velocity
% the incremental deformation of the element 
% Update max and min deformation history
% Element.prop.Dmax= max(Ue,Element.prop.Dmax); 
% Element.prop.Dmin= min(Ue,Element.prop.Dmin); 
du = Ue - Element.Uprev;
[Element.prop kt Fs] = MaterialState(du, Element.prop, Ve);

RestoringForce=zeros(6,1);
RestoringForce(Element.dof,1)   = Fs;
RestoringForce(Element.dof+3,1) = -Fs;

if( Structure.RigidLinkNo )
    RestoringForce = TR' * RestoringForce;  
end

% Update Uprev
Element.Uprev=Ue; 

