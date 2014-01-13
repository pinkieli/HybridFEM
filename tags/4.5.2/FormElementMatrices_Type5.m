%% Form global element stiffness matrix and mass matrix
%  ISTF = 1: Zero Matrices
%  ISTF = 0: Formed Matrices
function [ElementStiffnessMatrix,...
          ElementMassMatrix] = FormElementMatrices_Type5(Element, ISTF)
      
% Preset arrays
ElementStiffnessMatrix = zeros(6,6);
ElementMassMatrix = zeros(6,6);
    
if ISTF == 0    
    ElementStiffnessMatrix(3,3) = Element.elstk;
    ElementStiffnessMatrix(3,6) = -Element.elstk;
    ElementStiffnessMatrix(6,6) = Element.elstk;
    ElementStiffnessMatrix(6,3) = ElementStiffnessMatrix(3,6);
else   
    ElementStiffnessMatrix = zeros(6);     
end