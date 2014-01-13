%% Form global element stiffness matrix and mass matrix of elastic
%  beam-column element along with shear deformations
%  ISTF = 1: Zero Matrices
%  ISTF = 0: Formed Matrices
function [ElementStiffnessMatrix,...
          ElementMassMatrix] = FormElementMatrices_Type1(Element, ISTF)
      
% Get distance between nodes
X21 = Element.Nodes(2).Xcoord - Element.Nodes(1).Xcoord;
Y21 = Element.Nodes(2).Ycoord - Element.Nodes(1).Ycoord;
EL = sqrt(X21 * X21 + Y21 * Y21);
EAL = Element.Material.E * Element.Area / EL;
EIZL = Element.Material.E * Element.Inertia / EL;
% Calculation of density (the uniform load is only taken into account)
RO = abs(Element.ElemDistLoad) / 9.81;
beta = RO * EL / 420.0;
SEP = zeros(6);
MEP = zeros(6);
MEP(1,1) = 140;
MEP(1,2) = 0;
MEP(1,3) = 0;
MEP(1,4) = 70;
MEP(1,5) = 0;
MEP(1,6) = 0;
MEP(2,2) = 156;
MEP(2,3) = 22 * EL;
MEP(2,4) = 0.0;
MEP(2,5) = 54;
MEP(2,6) = -13 * EL;
MEP(3,3) = 4 * EL * EL;
MEP(3,4) = 0.0;
MEP(3,5) = 13 * EL;
MEP(3,6) = -3.0 * EL * EL;
MEP(4,4) = 140;
MEP(4,5) = 0.0;
MEP(4,6) = 0.0;
MEP(5,5) = 156;
MEP(5,6) = -22 * EL;
MEP(6,6) = 4.0 * EL * EL;
for I = 1:6
   for J = I:6
      MEP(J, I) = MEP(I, J);
   end
end
MEP=beta*MEP;
SEP(1,1) = EAL;
SEP(1,4) = -EAL;
SEP(4,4) = EAL;
SEP(2,2) = 12 * EIZL / EL ^ 2;
SEP(2,3) = 6 * EIZL / EL;
SEP(2,5) = -SEP(2, 2);
SEP(2,6) = SEP(2, 3);
SEP(3,3) = 4 * EIZL;
SEP(3,5) = -6 * EIZL / EL;
SEP(3,6) = 2 * EIZL;
SEP(5,5) = 12 * EIZL / EL ^ 2;
SEP(5,6) = -6 * EIZL / EL;
SEP(6,6) = 4 * EIZL;
for I = 1:6
   for J = I:6
      SEP(J, I) = SEP(I, J);
   end
end

DCOS(1,1) = X21 / EL;
DCOS(1,2) = Y21 / EL;
DCOS(1,3) = 0;
DCOS(2,1) = -DCOS(1, 2);
DCOS(2,2) = DCOS(1, 1);
DCOS(2,3) = 0;
DCOS(3,1) = 0;
DCOS(3,2) = 0;
DCOS(3,3) = 1;
ALMBDA = zeros(6);       
for K = 1:2
   IK = 3 * (K - 1);
   for I = 1:3
      for J = 1:3
         ALMBDA(I + IK, J + IK) = DCOS(I, J);
      end
   end
end
if ISTF == 0
    SE = SEP * ALMBDA;
    ME = MEP * ALMBDA;
    SEP = SE;
    MEP = ME;
    ElementStiffnessMatrix = ALMBDA' * SEP;
    ElementMassMatrix = ALMBDA' * MEP;
else
    ElementStiffnessMatrix = zeros(6);  
    ElementMassMatrix = zeros(6);  
end