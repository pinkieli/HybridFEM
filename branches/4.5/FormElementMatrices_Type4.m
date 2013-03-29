%% Form global element stiffness matrix and mass matrix of lean on column
%  ISTF = 1: Zero Matrices
%  ISTF = 0: Formed Matrices
function [ElementStiffnessMatrix,...
          ElementMassMatrix] = FormElementMatrices_Type4(Element, ISTF)
  
% dummy column element must share the nodes with lean on column
      
      
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

% construct geometric stiffness matrix and mass matrix in local coordinate
kele = zeros(6);
mass = zeros(6);
kele(2,2)=Element.K22;
kele(2,5)=Element.K25;
kele(5,2)=Element.K52;
kele(5,5)=Element.K55;
mass(:,:)=Element.m;
mass(Element.massDof,Element.massDof)=Element.m1; % mass on i node
mass(Element.massDof+3,Element.massDof+3)=Element.m2; % mass on j node

% 6x6 transformation matrix from local to global coordiante
ElementStiffnessMatrix = apq' * kele * apq;

% mass matrix in global coordinate
ElementMassMatrix =  mass ;

if (ISTF == 1)
    ElementStiffnessMatrix = zeros(6);  
    ElementMassMatrix = zeros(6);  
end