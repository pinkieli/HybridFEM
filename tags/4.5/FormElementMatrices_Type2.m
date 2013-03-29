%% Forms global element stiffness matrix and mass matrix (zero) of experimental link element...
%  assumes that the element is oriented in the X global direction and
%  connects two nodes with the same coordinates or different coordinates
%  ISTF = 1: Zero Matrices
%  ISTF = 0: Formed Matrices
function [ElementStiffnessMatrix,...
          ElementMassMatrix] = FormElementMatrices_Type2(Element, ISTF)
      
% Preset arrays
ElementStiffnessMatrix = zeros(6,6);
ElementMassMatrix = zeros(6,6);
    
if ISTF == 0    
    ElementStiffnessMatrix(1,1) = Element.Stiffness;
    ElementStiffnessMatrix(1,4) = -Element.Stiffness;
    ElementStiffnessMatrix(4,4) = Element.Stiffness;
    ElementStiffnessMatrix(4,1) = ElementStiffnessMatrix(1,4);
end