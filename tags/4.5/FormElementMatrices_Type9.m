%% Form global element stiffness matrix and mass matrix
%  ISTF = 1: Zero Matrices
%  ISTF = 0: Formed Matrices
function [ElementStiffnessMatrix,...
          ElementMassMatrix] = FormElementMatrices_Type9(Element, ISTF)
    
% don't know what ISTF is !!! just use it as Theodore did it for EType3
if ISTF == 0
    
    % element geometry       
    H = Element.Height; 
    L = Element.Width;
    
    dx = abs( Element.xyk(1,1)-  Element.xyi(1,1));
    dy = abs( Element.xyk(2,1)-  Element.xyi(2,1));
    
    % fill transformation coefficients, global coordinate to local coordinate
    avq = zeros(9,12);
    
    avq(1,1) =  -dx/L;
	avq(1,2) =  -dy/L;
    avq(1,7) =  dx/L;
	avq(1,8) =  dy/L;

    avq(2,4) =   dy/L;
	avq(2,5) =  -dx/L;
    avq(2,10)=  -dy/L;
	avq(2,11)=   dx/L;

    avq(3,1) =   H/L*dy/L;
	avq(3,2) =  -H/L*dx/L;
    avq(3,4) =  -dx/L;
	avq(3,5) =  -dy/L;    
    avq(3,7) =  -H/L*dy/L;
	avq(3,8) =   H/L*dx/L;
    avq(3,10)=   dx/L;
	avq(3,11)=   dy/L;   
    
    avq(4,3) =  -1.;
	avq(4,9) =   1.;
    avq(5,6) =  -1.;
	avq(5,12)=   1.;    

    avq(6,3) =   1.0;
	avq(6,4) =  -dx/L*2.0/H;
    avq(6,5) =  -dy/L*2.0/H;
	avq(6,9) =   1.0;
    avq(6,10)=   dx/L*2.0/H;
	avq(6,11)=   dy/L*2.0/H; 
    
    avq(7,1) =  -dy/L*2.0/L;
	avq(7,2) =  dx/L*2.0/L;
    avq(7,6) =  1.0;
	avq(7,7) =   dy/L*2.0/L;   
    avq(7,8) =  -dx/L*2.0/L;
	avq(7,12)=   1.0;
    
    avq(8,1) =   dx/L;
	avq(8,2) =   dy/L;
    avq(8,4) =  -dx/L;
	avq(8,5) =  -dy/L;    
    avq(8,7) =   dx/L;
	avq(8,8) =   dy/L;
    avq(8,10)=  -dx/L;
	avq(8,11)=  -dy/L;      
    
    avq(9,1) =   dy/L;
	avq(9,2) =  -dx/L;
    avq(9,4) =  -dy/L;
	avq(9,5) =   dx/L;    
    avq(9,7) =   dy/L;
	avq(9,8) =  -dx/L;
    avq(9,10)=  -dy/L;
	avq(9,11)=   dx/L;      
        
    % Integrate basic force and basic stiffness matrix
    kb = zeros(9);
   
    % basic stiffness matrix
    % calculate stiffness coefficient terms
    % make some stiffnees terms rigid  
    kb(1,1) = Element.kb1 ;
    kb(2,2) = Element.kb2 ;    
    kb(3,3) = Element.kb3 ;
    kb(4,4) = Element.kb4 ;    
    kb(5,5) = Element.kb5 ;
    kb(6,6) = Element.kb6 ;    
    kb(7,7) = Element.kb7 ;
    kb(8,8) = Element.kb8 ;    
    kb(9,9) = Element.kb9 ;    
    kb(3,3) = Element.ShearProp.Kt;  % Modified by Akbar
    % element stiffness matrix in global coordinate
    ElementStiffnessMatrix = avq' * kb * avq;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ElementMassMatrix=zeros(12);
    for i=1:12
        ElementMassMatrix(i,i) = Element.mass;
    end    
else 
    ElementStiffnessMatrix=zeros(12);
    ElementMassMatrix=zeros(12);
end

