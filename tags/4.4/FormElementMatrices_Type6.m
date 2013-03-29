%% Form global element stiffness matrix and mass matrix
%  ISTF = 1: Zero Matrices
%  ISTF = 0: Formed Matrices
function [ElementStiffnessMatrix,...
          ElementMassMatrix] = FormElementMatrices_Type6(Element, ISTF)
    
% element geometry       
% element deformation % geometry restricted to 1,2 plane (x,y)
dx = Element.xyj(1,1) - Element.xyi(1,1); 
dy = Element.xyj(2,1) - Element.xyi(2,1);
L = sqrt(dx*dx+dy*dy);
dx = dx / L;
dy = dy / L;

% Form transformation matrix from element in global coordinate to basic system 
avq  = [ -dy/L  dx/L  1  dy/L -dx/L  0;
        -dy/L  dx/L  0  dy/L -dx/L  1;
        -dx   -dy    0  dx    dy    0];
% Incremental basic deformations
dv = avq * [0 0 0 0 0 0]';

% Interpolation functions (B matrix) evaluated at integration points
x  = ( Element.xi + 1 ) / 2;
b1 = (6*x-4)/L;
b2 = (6*x-2)/L;
b3 = ones(size(x))/L;

% Integrate basic force and basic stiffness matrix
kb = zeros(3);
s  = zeros(3,1);

for i=1:length(Element.sections)
    % Form B matrix and compute section deformation
    B = [ b1(i) b2(i) 0; 0 0 b3(i) ];
    dvs = B * dv;  

    % Section state determination and update section state
    [ss ks sec] = SectionState( Element.sections(i), dvs );

    % Sum section term
    kb = kb + B'*ks*B * ( Element.wi(i) * L/2 );
    s  = s  + B'*ss   * ( Element.wi(i) * L/2 );
end

% element stiffness matrix in global coordinate
ElementStiffnessMatrix = avq' * kb * avq;


% Form consistent mass matrix, copied from FormElementMatrices_Type5
% Calculation of density (the uniform load is only taken into account)
rho = abs(Element.ElemDistLoad) / 9.81;
beta = rho * L / 420.0;

% consistent mass matrix 
m = zeros(6);
m(1,1) = 140.;
m(1,2) = 0;
m(1,3) = 0;
m(1,4) = 70.;
m(1,5) = 0;
m(1,6) = 0;
m(2,2) = 156.;
m(2,3) = 22. * L;
m(2,4) = 0.0;
m(2,5) = 54.;
m(2,6) = -13. * L;
m(3,3) = 4. * L * L;
m(3,4) = 0.0;
m(3,5) = 13. * L;
m(3,6) = -3.0 * L * L;
m(4,4) = 140.;
m(4,5) = 0.0;
m(4,6) = 0.0;
m(5,5) = 156.;
m(5,6) = -22. * L;
m(6,6) = 4.0 * L * L;
for I = 1:6
    for J = I:6
        m(J, I) = m(I, J);
    end
end
m=beta*m;

% 6x6 transformation matrix from local to global coordiante
apq  = zeros(6);
apq(1,1) = dx; apq(2,2) = dx; apq(4,4) = dx; apq(5,5)= dx;
apq(1,2) = dy; apq(2,1) =-dy; apq(4,5) = dy; apq(5,4)=-dy;
apq(3,3) =1.0; apq(6,6) =1.0;
% element mass matrix in global coordinate
ElementMassMatrix = apq' * m * apq;

if ISTF == 1
    ElementStiffnessMatrix=zeros(6);
    ElementMassMatrix=zeros(6);
end
