%% Create Element Type 6 Displacement Based Beam-Column Element
function vars = CreateElement_Type6(ID, Type, Node1, Node2, Materials,...
    DampStiffFac, DampMassFac, Section, nIP, ElemDistLoad)

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
vars.L = L;
if (L < eps)
	error(['Element number ' ID ' has zero length.']);
end
    
% assign sections to the element
for i=1:nIP
    vars.sections(i) = Section;   
end
dx = dx / L;
dy = dy / L;

% Form transformation matrix from element in global coordinate to basic system 
vars.avq  = [ -dy/L  dx/L  1  dy/L -dx/L  0;
        -dy/L  dx/L  0  dy/L -dx/L  1;
        -dx   -dy    0  dx    dy    0];

% section locations and integration variables
[vars.xi, vars.wi] = GaussLobatto( nIP );

% Interpolation functions (B matrix) evaluated at integration points
x  = ( vars.xi + 1 ) / 2;
vars.b1 = (6*x-4)/L;
vars.b2 = (6*x-2)/L;
vars.b3 = ones(size(x))/L;

% Set handles to element matrices forming function and restoring force
vars.FormElementMatrices = @FormElementMatrices_Type6;
vars.GetRestoringForce = @GetRestoringForce_Type6;

% Keep track of element displacement vector in global coordinates
vars.Uprev=zeros(6,1);


% function [xi, wi] =  Gauss(nIP)
% % locations (xi) and weights (wi) for the number of integration points
% 
% switch nIP
% 
%     case 2
%         xi= [ -.57735026918963,.57735026918963]';
%     	wi= [ 1.,1.]';
% 
%   	 case 3
%     	xi= [-0.77459666924148,0.0,0.77459666924148]';
%     	wi= [.55555555556,.88888888889,.55555555556]';
%     
%     case 4
%         xi= [-.8611363116,-.3399810436,.3399810436,.8611363116]';
%     	wi= [.3478548451,.6521451549,.6521451549,.3478548451]';
%  
% 	case 5
%         xi= [-.9061798459,-.5384693101,0.0,.5384693101,.9061798459]';
%     	wi= [.236926885,.4786286705,.5688888889,.4786286705,.236926885]';
% 
% 	case 6
%         xi= [-.9324695142,-.6612093865,-.2386191861, ...
%        			 .2386191861, .6612093865, .9324695142]';
%         wi= [.1713244924,.3607615730,.4679139346, ...
%       			 .4679139346,.3607615730,.1713244924]';
%     case 7
%         xi= [-.9491079123,-.7415311856,-.4058451514,0., ...
%        			 .4058451514, .7415311856, .9491079123]';
%         wi= [.1294849662,.2797053915,.3818300505,.4179591837, ...
%       			 .3818300505,.2797053915,.1294849662]';
%     otherwise
%         error('Invalid order for Gauss integration.')
% 	
% end

function [xi, wi] =  GaussLobatto(nIP)
% this lobatto integration scheme is used in OpenSEES comparable element
% by default. So we are using the same scheme for comparison purpose

% locations (xi) and weights (wi) for the number of integration points
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