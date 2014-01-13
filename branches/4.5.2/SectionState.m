%% State determination for fiber section
function [ss, ks, sec] = SectionState( sec, dvs )
% For section sec undergoing a change
% in section deformation dvs=[curvature strain], compute
% section forces SS=[moment axial] and the 2x2 tangents
% section stiffness KS.

ks = zeros(2,2);
ss = zeros(2,1);  
numfib = length(sec.fibers);

% Loop over fibers
for i=1:numfib
	
	yi   = sec.fibers(i).yi;
	area = sec.fibers(i).A;
	
	% Change in strain
	deps = yi * dvs(1:2);
	
    [sec.fibers(i).prop Et sc] = MaterialState(deps, sec.fibers(i).prop, 0 );
    %sc = sec.fibers(i).prop.sc;

	% Integration -- add contribution of fiber 
    % for 3D --- following codes need to be modified
	ks = ks + yi'*yi * (area*Et ); % may not need to do it
	ss = ss + yi'    * (area*sc);

end


