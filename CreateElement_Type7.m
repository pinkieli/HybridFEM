%% Create Element Type 7 
% flexiblity based beam-column element 
function vars = CreateElement_Type7(ID, Type, Node1, Node2, Materials,...
    DampStiffFac, DampMassFac, Section, nIP, ElemDistLoad, maxIter, EAT,ERT,CO)

% Section geometry and properties read from input file
% Section(1)  = Section ID
% Section(2)  = Section Type
% Section(3)  = Section height, h
% Section(4)  = Flange width, bf
% Section(5)  = Flange thickness, tf
% Section(6)  = Web width, tw
% Section(7)  = No. fibers per flange, nf
% Section(8)  = No. fibers for web
% Section(9)  = material type 
% if Section(9) == 1 , bilinear material 
%    Section(10) = E, Section(11) = sigmaY, Section(12) = second slop ratio
% if Section(9) == 3 , concrete material  
% CO: If you carry over the unb residual element disp= 1, if not CO=0

% Set element nodes
vars.ID = ID;
vars.Type=Type;
vars.Nodes = [Node1 Node2];
vars.DampStiffFac=DampStiffFac;
vars.DampMassFac=DampMassFac;
vars.ElemDistLoad = ElemDistLoad;

% Variables which are initialized from the input
vars.xyj(1,1) =  Node2.Xcoord;    
vars.xyj(2,1) =  Node2.Ycoord;
vars.xyi(1,1) =  Node1.Xcoord;
vars.xyi(2,1) =  Node1.Ycoord;

% calculate the element length
dx = vars.xyj(1,1) - vars.xyi(1,1);
dy = vars.xyj(2,1) - vars.xyi(2,1);
L = sqrt(dx*dx+dy*dy);   
if (L < eps)
	error(['Element number ' ID ' has zero length.']);
end
    
% assign sections to the element
vars.nIP = nIP;
for i=1:nIP
    vars.sections(i) = Section;   
end

% convert and copy the sections to a matrix form for use with the simulink element version
%vars.SectionMatrix = CopySectionsToMatrix(vars.sections, Materials);

% section locations and integration variables
[vars.xi, vars.wi] = Lobatto( nIP );

% criteria to satisfy compatibility during iterative procedure
vars.maxIter = maxIter;
vars.NumIter=1;
vars.EAT = EAT;
vars.ERT1 = ERT;
vars.EC=0;
vars.CO = CO;
vars.NUnbforce=0;
% Set handles to element matrices forming function and restoring force
vars.FormElementMatrices = @FormElementMatrices_Type7;
vars.GetRestoringForce = @GetRestoringForce_Type7;
    
% Keep track of element displacement vector in global coordinates
vars.Uprev=zeros(6,1);
vars.resdq=zeros(6,1);
vars.si = zeros(3,1);


function [xi, wi] =  Lobatto(nIP)
% locations (xi) and weights (wi) for the number of integration points
xi    = zeros(nIP,1);
wi    = zeros(nIP,1);

switch nIP
    case 2
        %xi= [-1., 1.]';  %Gauss-Labatto
        xi  = [-0.577350269189626 0.577350269189626]';  %Gauss-Legendre
        wi = [ 1. 1.]';
        %errordlg('Gauss-Legendre is used for integration points- RC-MRF Analysis ')
    case 3
        %xi= [-1.,0., 1.]'; %Gauss-Labatto
    	%wi= [1/3, 1/3, 1/3]'; %Gauss-Labatto
        xi= [-0.7745966692,0., 0.7745966692]'; %Gauss-Legendre
    	wi= [0.5555555556, 0.8888888888, 0.5555555556]'; %Gauss-Legendre
	case 4
        %xi= [-1.,-0.44721360, 0.44721360, 1.]'; %Gauss-Labatto
    	%wi= [1/6, 5/6, 5/6, 1/6]'; %Gauss-Labatto
        xi= [-0.8611363116,-0.3399810436, 0.3399810436, 0.8611363116]'; %Gauss-Legendre
    	wi= [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]'; %Gauss-Legendre
	case 5
        %xi= [-1.,-0.65465367, 0., 0.65465367, 1.]';  %Gauss-Labatto
    	%wi= [0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1]';  %Gauss-Labatto
        xi= [-0.9061798459,-0.5384693101, 0., 0.5384693101, 0.9061798459]';  %Gauss-Legendre
    	wi= [0.2369268851, 0.4786286705, 0.5688888888, 0.4786286705, 0.2369268851]';  %Gauss-Legendre
    case 6
        xi= [-1.,-0.7650553239, -0.2852315164, 0.2852315164,...
            0.7650553239, 1.]';
    	wi= [0.06666666667, 0.3784749562, 0.5548583770, 0.5548583770,...
            0.3784749562, 0.06666666667]';
    case 7
        xi= [-1.,-0.8302238962, -0.4688487934, 0.0, 0.4688487934,...
            0.8302238962, 1.]';
    	wi= [0.04761904762, 0.2768260473, 0.4317453812, 0.4876190476,...
            0.4317453812, 0.2768260473, 0.04761904762]';
    otherwise
        error('Invalid order for Lobatto integration.')
end
	