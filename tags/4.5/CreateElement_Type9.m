%% Create Element Type 9 
% Flexiblity based beam-column element 
function vars = CreateElement_Type9(ID, Type, Node1, Node2, Node3, Node4,...
    DampStiffFac, DampMassFac, mat, Icl, Acl, bfcl, twcl, tfcl, tdp, tcnt,...
    shearMat, mass)
            
% Set element nodes
vars.ID = ID;
vars.Type=Type;
vars.Nodes = [Node1 Node2 Node3 Node4];
vars.DampStiffFac=DampStiffFac;
vars.DampMassFac=DampMassFac;

% Variables which are initialized from the input
vars.xyi(1,1) =  Node1.Xcoord;
vars.xyi(2,1) =  Node1.Ycoord;
vars.xyj(1,1) =  Node2.Xcoord;    
vars.xyj(2,1) =  Node2.Ycoord;
vars.xyk(1,1) =  Node3.Xcoord;    
vars.xyk(2,1) =  Node3.Ycoord;
vars.xyl(1,1) =  Node4.Xcoord;
vars.xyl(2,1) =  Node4.Ycoord;

vars.Height = abs(Node2.Ycoord - Node4.Ycoord);
vars.Width  = abs(Node1.Xcoord - Node3.Xcoord);

% Check element dimension
if vars.Height <= 1e-12 || vars.Width <= 1e-12
    error(['Element number ' ID ' has either zero length or width']);
end

% transform displacement to deformation 
vars.dx = abs( vars.xyk(1,1)-  vars.xyi(1,1));
vars.dy = abs( vars.xyk(2,1)-  vars.xyi(2,1));

ndof=12;  
avq = zeros(9,ndof);

avq(1,1) =  -vars.dx/vars.Width;
avq(1,2) =  -vars.dy/vars.Width;
avq(1,7) =   vars.dx/vars.Width;
avq(1,8) =   vars.dy/vars.Width;

avq(2,4) =   vars.dy/vars.Width;
avq(2,5) =  -vars.dx/vars.Width;
avq(2,10)=  -vars.dy/vars.Width;
avq(2,11)=   vars.dx/vars.Width;

avq(3,1) =   vars.Height/vars.Width*vars.dy/vars.Width;
avq(3,2) =  -vars.Height/vars.Width*vars.dx/vars.Width;
avq(3,4) =  -vars.dx/vars.Width;
avq(3,5) =  -vars.dy/vars.Width;    
avq(3,7) =  -vars.Height/vars.Width*vars.dy/vars.Width;
avq(3,8) =   vars.Height/vars.Width*vars.dx/vars.Width;
avq(3,10)=   vars.dx/vars.Width;
avq(3,11)=   vars.dy/vars.Width;   

avq(4,3) =  -1.;
avq(4,9) =   1.;
avq(5,6) =  -1.;
avq(5,12)=   1.;    

avq(6,3) =   1.0;
avq(6,4) =  -vars.dx/vars.Width*2.0/vars.Height;
avq(6,5) =  -vars.dy/vars.Width*2.0/vars.Height;
avq(6,9) =   1.0;
avq(6,10)=   vars.dx/vars.Width*2.0/vars.Height;
avq(6,11)=   vars.dy/vars.Width*2.0/vars.Height; 

avq(7,1) =  -vars.dy/vars.Width*2.0/vars.Width;
avq(7,2) =  vars.dx/vars.Width*2.0/vars.Width;
avq(7,6) =  1.0;
avq(7,7) =   vars.dy/vars.Width*2.0/vars.Width;   
avq(7,8) =  -vars.dx/vars.Width*2.0/vars.Width;
avq(7,12)=   1.0;

avq(8,1) =   vars.dx/vars.Width;
avq(8,2) =   vars.dy/vars.Width;
avq(8,4) =  -vars.dx/vars.Width;
avq(8,5) =  -vars.dy/vars.Width;    
avq(8,7) =   vars.dx/vars.Width;
avq(8,8) =   vars.dy/vars.Width;
avq(8,10)=  -vars.dx/vars.Width;
avq(8,11)=  -vars.dy/vars.Width;      

avq(9,1) =   vars.dy/vars.Width;
avq(9,2) =  -vars.dx/vars.Width;
avq(9,4) =  -vars.dy/vars.Width;
avq(9,5) =   vars.dx/vars.Width;    
avq(9,7) =   vars.dy/vars.Width;
avq(9,8) =  -vars.dx/vars.Width;
avq(9,10)=  -vars.dy/vars.Width;
avq(9,11)=   vars.dx/vars.Width; 

vars.avq = avq;
    
% Panel zone material and section properties
E = mat.E;
G = E/(2.*(1.+0.3)); 

d1 = vars.Width ;
d2 = vars.Height + 2.0*tcnt;

%% Element stiffness coefficients
%  Axial stiffness
SecA1    = d2*(twcl+tdp) + 2.*tcnt*(bfcl-twcl) ;
vars.kb1 = E*SecA1/d1;
SecA2    = Acl + (d1 - 2.*tfcl)*tdp;
vars.kb2 = E*SecA2/d2;
%  Bending stiffness
I11 = d2*d2*d2*twcl/12.;
I12 = 2.*tcnt*tcnt*tcnt*(bfcl-twcl)/12.;
I13 = 2.*(bfcl-twcl)*tcnt*(d2-tcnt)*(d2-tcnt)/4;
I14 = d2*d2*d2*tdp/12.;
vars.kb4 = E*(I11+I12+I13+I14)/d1;
I21 =  (d1-2.*tfcl)*(d1-2.*tfcl)*(d1-2.*tfcl)*tdp/12.;
vars.kb5 = E*(I21+Icl)/d2;
%  Asym Bending stiffness
vars.kb6 = 3.*E*(I11+I12+I13+I14)/d1;
vars.kb7 = 3.*E*(I21+Icl)/d2;
%  Asym Shear Stiffness
vars.kb8 = 6.*E*G*SecA1*d1/(3.*E*d2*d2 +2.*G*d1*d1);    
vars.kb9 = 6.*E*G*SecA2*d2/(3.*E*d1*d1 +2.*G*d2*d2);    
vars.kb3 = G*(twcl+tdp)*d1/d2;

%% make some deformation modes rigid
vars.kb1 = vars.kb1*100.;
vars.kb2 = vars.kb2*100.;
vars.kb4 = vars.kb4*3.;
vars.kb5 = vars.kb5*3.;
vars.kb6 = vars.kb6*100.;
vars.kb7 = vars.kb7*100.;
vars.kb8 = vars.kb8*100.;    
vars.kb9 = vars.kb9*100.;    

% Element shear properties
vars.ShearProp = shearMat;   

% Mass associated with panel zone, typically small mass 
vars.mass = mass;

% Set handles to element matrices forming function and restoring force
vars.FormElementMatrices = @FormElementMatrices_Type9;
vars.GetRestoringForce = @GetRestoringForce_Type9;

% Keep track of element displacement vector in global coordinates
vars.Uprev=zeros(12,1);
vars.si = zeros(9,1);
