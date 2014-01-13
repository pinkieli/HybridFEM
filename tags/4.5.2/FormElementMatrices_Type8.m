%% Form global element stiffness matrix and mass matrix
%  ISTF = 1: Zero Matrices
%  ISTF = 0: Formed Matrices
function [ElementStiffnessMatrix,...
          ElementMassMatrix] = FormElementMatrices_Type8(Element, ISTF)
      
% Preset arrays
ElementStiffnessMatrix = zeros(6,6);
ElementMassMatrix = zeros(6,6);
    
if ISTF == 0    
    
    % element stiffness ---- later change material input parameters to have
    % consistent notation so we dont need following switch statement
    switch Element.prop.Type
        case 1
            K = Element.prop.E; % for elastic material Material
        case 2
            K = Element.prop.e; % for Bilinear material
        case 3
            K = Element.prop.E1p; % for Hysteretic material
        case 4
            K = Element.prop.k; % for Bouc-Wen Material
        case 5
            K = Element.prop.Kt; % for Trilinear material
        case 6
            K = Element.prop.Kip; % for SDegrading material
        otherwise
            errordlg('The assigned material is not available', 'InputError');
    end
    
    ElementStiffnessMatrix(Element.dof , Element.dof ) = K;
    ElementStiffnessMatrix(Element.dof , Element.dof+3)= -K;    
    ElementStiffnessMatrix(Element.dof+3,Element.dof+3) = K;
    ElementStiffnessMatrix(Element.dof+3,Element.dof )= -K;    
    
else
    ElementStiffnessMatrix=zeros(6);
    ElementMassMatrix=zeros(6);
end