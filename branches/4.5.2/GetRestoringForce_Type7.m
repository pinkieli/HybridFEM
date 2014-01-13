%% Get restoring force for the element
%  returns (6x1) restoring force vector of element in global coordinates
function [Element,RestoringForce] = GetRestoringForce_Type7(Element, Structure, Integrator)

% Get the displacement Ue (6x1) vector of the element in global coordinates
Ue = zeros(6,1);
dof=zeros(6,1);
dof(1,1)=Element.Nodes(1).UX;
dof(2,1)=Element.Nodes(1).UY;
dof(3,1)=Element.Nodes(1).THETA;
dof(4,1)=Element.Nodes(2).UX;
dof(5,1)=Element.Nodes(2).UY;
dof(6,1)=Element.Nodes(2).THETA;
for i=1:6
    if dof(i,1) ~= -1
        Ue(i,1)=Integrator.Displacement(dof(i,1),1);
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
    Ue = TR * Ue;   
end

ndof=6;  % The element has 6 dofs

if Element.CO==1
    du=(Ue+Element.resdq)-Element.Uprev;
else
    du=Ue -Element.Uprev; % The incremental displacement vector du of the element in global coordinates 
end

Element.Uprev=Ue; %Update Uprev of element for use in next state determination call

% element defomration % geometry restricted to 1,2 plane (x,y)
dx = Element.xyj(1,1) - Element.xyi(1,1); 
dy = Element.xyj(2,1) - Element.xyi(2,1);
L = sqrt(dx*dx+dy*dy);
dx = dx / L;
dy = dy / L;

% Form transformation matrix from element in global coordinate to basic system
avq  = [ -dy/L  dx/L  1.  dy/L -dx/L  0;
       -dy/L  dx/L  0  dy/L -dx/L  1.;
       -dx   -dy    0  dx    dy    0];

% Incremental basic deformations
dv = avq * du;

% Basic forces from previous load step
si = Element.si; %momentI, momentJ, axial

% Initialize ds for iteration
ds = zeros(3,1); %deltaMomentI, deltaMomentJ, deltaAxial

% Interpolation functions for equilibrium evaluated at integration points
x  = (  Element.xi + 1. ) / 2.;   
b1 = -( 1. - x );
b2 =        x;      
b3 = -1./L;
b4 = ones(size(x));

% initialize incremental section deformations
dvsecs = zeros(2,Element.nIP);

%vr = dv;
	
% Iterate to satisfy compatibility 
for k =1:Element.maxIter
  	
    % Compute residual and basic flexibility by summing 
    % section properties over element length
  	vr = dv;
	fb = zeros(3);
    si = Element.si; %momentI, momentJ, axial
    % loop over each section along the element length
    for i=1:Element.nIP
       
        % transformaion matrix 
        b = [ b1(i) b2(i) 0;
               0     0    b4(i)];
       %        b3    b3    0]; 
		% section force, moment, axial, and shear 
        ss = b * ( ds + si );
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Given section force, determine section deformation iteratively
        % Initialize vector of dvs for iteration   
         dvs = zeros(2,1);

        % Iterate on section equilibrium
  %      for j=1:Element.maxIter
	    
            % given section deformation, determine section force
            [jthss, ks] = SectionState(Element.sections(i), dvs);            
            % residual section force 
            sres     = ss - jthss;
            dvs = ks \ sres + dvs;    
            


  %      end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store deformation of each section
        dvsecs(:,i) = dvs;				
		% add section contribution to residual element deformation and flexibility
		vr = vr - b' * dvs   * ( Element.wi(i) * L/2. );
		fb = fb + b'*(ks\b) * ( Element.wi(i) * L/2. ); 
        
      
    end
 
    % check residual deformation 
%     TolFinal=min(Element.EAT,Element.ERT1*norm(vr));
%     if norm(vr) < TolFinal
%         Element.NumIter=k;
%         Element.EC=Element.EC+1;
% 		break;
%     end

% check Unbalanced force   
  UnbForce=fb\vr;
  NormUnb=((UnbForce(1)*L)^2+(UnbForce(2)*L)^2+(UnbForce(3))^2)^0.5;
  TolFinal=max(Element.EAT,Element.ERT1*NormUnb);
    if NormUnb < TolFinal
        Element.NumIter=k;
        Element.EC=Element.EC+1;
        Element.NUnbforce=NormUnb;
		break;
    end

    ds = fb \ vr + ds; 
 Element.NumIter=Element.maxIter;  
end

% update section state with converged dvs
for i=1:Element.nIP
	[ssc, ksc, Element.sections(i)] = SectionState(Element.sections(i),...
        dvsecs(:,i));    
end

% update basic forces for element
Element.si = ds + si;
if Element.CO==1
    Element.resdq = avq'*vr;    
end
% Restoring force and stiffness matrix
RestoringForce = avq' * Element.si;

if( Structure.RigidLinkNo )
    RestoringForce = TR' * RestoringForce;  
end