%% Get restoring force for the element
%  returns (6x1) restoring force vector of element in global coordinates
function [Element,RestoringForce] = GetRestoringForce_Type9(Element, Structure, Integrator)

%% Get the displacement Ue (12x1) vector of the element in global coordinates
ndof=12;  
Ue = zeros(ndof,1);
dof=zeros(ndof,1);
dof(1,1)=Element.Nodes(1).UX;
dof(2,1)=Element.Nodes(1).UY;
dof(3,1)=Element.Nodes(1).THETA;
dof(4,1)=Element.Nodes(2).UX;
dof(5,1)=Element.Nodes(2).UY;
dof(6,1)=Element.Nodes(2).THETA;
dof(7,1)=Element.Nodes(3).UX;
dof(8,1)=Element.Nodes(3).UY;
dof(9,1)=Element.Nodes(3).THETA;
dof(10,1)=Element.Nodes(4).UX;
dof(11,1)=Element.Nodes(4).UY;
dof(12,1)=Element.Nodes(4).THETA;

for i=1:ndof
    if dof(i,1) ~= -1
        Ue(i,1)=Integrator.Displacement(dof(i,1),1);
    end
end

if( Structure.RigidLinkNo )
    TR = eye(12,12); % Transformation matrix for rigid link
    for jjj=1:4
        idx=find(Structure.RigidLinkNodeID(:,2) == Element.Nodes(jjj).ID); % find if a slave node exists in the element
        if ( idx )
            TL=CalculateRigidLinkTransformation(Structure.RigidLinkMaster{idx},Structure.RigidLinkSlave{idx} );
            TR(3*(jjj-1)+1:3*jjj,3*(jjj-1)+1:3*jjj) = TL;
        end
    end     
    Ue = TR * Ue;     
end

% The incremental displacement vector du of the element in global coordinates
% The element has 12 dofs
% Update Uprev of element for use in next state determination call
du=Ue-Element.Uprev;  
Element.Uprev=Ue; 

% transform global node displacement to basic deformation
v = Element.avq*Ue;

% compute element force
s = zeros(9,1); 

s(1) = Element.kb1*v(1);
s(2) = Element.kb2*v(2);
s(4) = Element.kb4*v(4);
s(5) = Element.kb5*v(5);
s(6) = Element.kb6*v(6);
s(7) = Element.kb7*v(7);
s(8) = Element.kb8*v(8);
s(9) = Element.kb9*v(9);

dv3 = Element.dx/Element.Width*du(10)-Element.dy/Element.Width*du(5)-Element.dx/Element.Width*du(4)+Element.dy/Element.Width*du(11)-Element.Height/Element.Width*Element.dx/Element.Width*du(2)...
    +Element.Height/Element.Width*Element.dy/Element.Width*du(1)+Element.Height/Element.Width*Element.dx/Element.Width*du(8)-Element.Height/Element.Width*Element.dy/Element.Width*du(7);

% compute shear deformation mode
[Element.ShearProp , ~, s(3)]= MaterialState(dv3, Element.ShearProp);

% Restoring force and stiffness matrix
RestoringForce = Element.avq' * s;

if( Structure.RigidLinkNo )
    RestoringForce = TR' * RestoringForce;  
end

% function [pro] = trilinear(dx, pro)
% 
% % calculate current deformation 
%     x = pro.xprev + dx; 
% 
% if pro.YieldCode == 0
% 
% 	pro.Fs1 = pro.Fs1+pro.K1*dx;
%     pro.Fs2 = pro.Fs2+pro.K2*dx;
% 	% if it is beyond elastic range, limit to yield force
% 	% and set m_nYieldCode to 1
% 	if pro.Fs1 > pro.Fy1
% 		
%         pro.Fs1 = pro.Fy1 ;
% 		pro.YieldCode = 1;
%         
%         if(pro.Fs2> pro.Fy2)
% 			pro.Fs2 = pro.Fy2 ;
% 		    pro.YieldCode = 2;
%         end
%     elseif pro.Fs1 < -pro.Fy1
% 		
%         pro.Fs1 = -pro.Fy1;
% 		pro.YieldCode = 1 ;
%        
% 		if pro.Fs2 < -pro.Fy2
% 			pro.Fs2 = - pro.Fy2 ;
% 		    pro.YieldCode = 2;
%         end
%     end
% % 1st post-yielding region    
% elseif pro.YieldCode == 1
%     
%     % assume spring2 remain elastic
%     pro.Fs2= pro.Fs2+pro.K2*dx;
%     
%     % keep loading
%     if pro.Fs1*dx > 0.
%         if pro.Fs2 > pro.Fy2
%             pro.Fs2 = pro.Fy2;
%             pro.YieldCode=2;            
%         elseif pro.Fs2 < -pro.Fy2
%             pro.Fs2 = -pro.Fy2;
%             pro.YieldCode=2;
%         end
%         % unloading
%     elseif pro.Fs1*dx < 0.
%         pro.Fs1= pro.Fs1 + pro.K1*dx;
%         pro.YieldCode = 0;
%         
%         if pro.Fs1 < -pro.Fy1 
%             pro.Fs1 = -pro.Fy1;
%             pro.YieldCode = 1;
%             if pro.Fs2 < -pro.Fy2
%                 pro.Fs2 = -pro.Fy2;
%                 pro.YieldCode = 2;
%             end
%         elseif pro.Fs1 > pro.Fy1
%             pro.Fs1 = pro.Fy1;
%             pro.YieldCode = 1;
%             if pro.Fs2 > pro.Fy2
%                 pro.Fs2 = pro.Fy2;
%                 pro.YieldCode = 2;
%             end
%         end
%     end
% elseif pro.YieldCode ==2
%     % unloading
% 	if (pro.Fs1*dx) < 0.
% 		pro.Fs1 = pro.Fs1 + pro.K1*dx;
% 		pro.Fs2 = pro.Fs2 + pro.K2*dx;
% 		pro.YieldCode = 0;
% 		
%         if pro.Fs1 < -pro.Fy1
%             pro.Fs1 = -pro.Fy1;
% 			pro.YieldCode = 1;
%             if pro.Fs2 < -pro.Fy2
% 				pro.Fs2 = -pro.Fy2;
% 				pro.YieldCode = 2;
%             end            
%         elseif pro.Fs1 > pro.Fy1
% 			pro.Fs1 = pro.Fy1;
% 			pro.YieldCode = 1;
%             if pro.Fs2> pro.Fy2
% 				pro.Fs2 = pro.Fy2;
% 				pro.YieldCode = 2;
% 			end
%         end
%     end
% end
% 
% pro.Fs = pro.Fs1 + pro.Fs2 + pro.K3*x;
% 
% % updata current state for next time step
% if pro.YieldCode == 0
%     pro.Kt = pro.K1 +pro.K2 +pro.K3;
% elseif pro.YieldCode == 1
%     pro.Kt = pro.K2 +pro.K3;    
% else
%     pro.Kt = pro.K3;
% end
% pro.xprev = x;